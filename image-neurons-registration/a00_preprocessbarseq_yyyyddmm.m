%% Preprocess images
clc

p=gcp('nocreate');
if isempty(p)
    parpool("local",12);
elseif p.NumWorkers<12
    parpool("local",12);
end
%%
organize_geneseq;
% run n2v in python
clc
% copyfile processed processed-copy
cmdout=n2v_geneonly();

% make a copy, just in case.

copyfile processed processed-n2vcopy
cd processed 

% Register all cycles. 
% This sections registers all sequencing cycles individuall, and register 40x BC to 40x genes
clc

%% Load channel profile for 20x
load chprofile20x-50-30-20-40-20220218.mat
load chshift20x-20220218
chshift20x=chshift20x(1:4,:);%in case this contains 5 channels.
chshift20x=chshift20x+[1,-0.7;1,0;0,0;0,0]; %manually fixed based on current data.

fname='n2vgene';
ball_radius=6;
rgb_intensity=0.3;
tic
register_seq_images(fname,chprofile20x,ball_radius,chshift20x,rgb_intensity)
t=toc;
%%
%remake_rgb_img(fname,rgb_intensity)

%94 hours 36minutes
%% Make codebook for genes
make_codebook('codebookM1all.mat')
%
% Process nuclear stain and sequential rounds. Before running this, start bardensr on GPU for basecalling, then run this at the same time.
% organize_hyb_files('hyb01')

% register hyb to seq first cycle
% load/make chprofile for sequential cycles
%load chprofilehyb
clc
chprofilehyb=[1 0 0 0;0.52 1 0 0;0 0 1 0;0 0 0 1]; % measured on this dataset
chshifthyb=[0 -3;0 0;0 0;-1 4]; % measured on this dataset
save('chprofilehyb.mat','chprofilehyb','chshifthyb');
%% 

% align sequential image to first sequencing cycle using the sequencing signals.
ch=3; %channel used to align
nuclearch=5; %nuclear channel
radius=100; %radius for bgn subtraction of all channels.
reg_cycle=1; % register to 1st cycle
chradius=30; % bgn radius for seq rolony cycle
rgb_intensity=0.3;

register_hyb_images(ch,radius, nuclearch,reg_cycle,chradius,chprofilehyb,chshifthyb, rgb_intensity);


%% Stitch images from first cycle
ch_count=5;
reg_ch=4;
overlap=0.23;
x=1;
trim=0.02;
rescale_factor=0.5; % rescaling images by this

stitch_10x_images(ch_count,reg_ch,overlap,x,trim,rescale_factor);

%save temp
%%
% Here check the images in "stitched40xto10x" to see if the two images match.        
check_stitching('40xto10x.mat')

% stitch the transformed 40x image and overlay with 10x images for visual inspection. This produces 5x smaller images than the original 10x to ease loading, but should still be enough since the registration only need to be roughly correct.

% The next section basecall genes.


% stitch the transformed 40x image and overlay with 10x images for visual inspection. This produces 5x smaller images than the original 10x to ease loading, but should still be enough since the registration only need to be roughly correct.
% check registration for all cycles
clc
intensity_scaling=3;
checkregistration('geneseq','40xto10x.mat',intensity_scaling);
% checkregistration('nuclear','40xto10x.mat',intensity_scaling);


% Here check the images in "chekregistration" to see if all the images match within each pos.        
% The next section basecall genes.



%% ORganize barcodes. This can be done either here or in a01 depending on whether data processing is concurrent with data collection

%move files

clc
cd ..
organize_bcseq()
%%
clc
cmdout=n2v_processing('n2vprocessing_bc.py');


%BCseq01 to geneseq01
%
%% 


cd processed

BC_refch=5;
gene_refch=5;
BC_name='n2vbcseq';
gene_name='n2vgeneseq';

alignBC2gene(BC_refch,gene_refch,BC_name,gene_name);

%% Register BC 40x

clc

% Load channel profile for 20x
load chprofile20x-50-30-20-40-20220218.mat
load chshift20x-20220218
chshift20x=chshift20x(1:4,:);%in case this contains 5 channels.
chshift20x=chshift20x+[1,-0.7;1,0;0,0;0,0]; %manually fixed based on current data.

fname='regn2vbc';
ball_radius=100;
rgb_intensity=0.3;


register_seq_images(fname,chprofile20x,ball_radius,chshift20x,rgb_intensity)
%
intensity_scaling = 3;

checkregistration('bcseq','40xto10x.mat',intensity_scaling);

% Stitch BC image, this is optional, to produce a first-cycle BC
%stitched image
niestitchgenessingleplane(5,4,0.23,1,0.02,'regn2vbc');
%%
checkregistration('bcseq','40xto10x.mat',intensity_scaling);
%%

function cmdout=n2v_processing(fname)
    if ~exist('fname','var')
        fname='n2vprocessing.py';
    end

    %This needs to be changed on new systems
    py_root_n2v = fileparts("C:\Users\mara.rue\Anaconda3\envs\n2v_mara\python.exe");
    ENV = getenv('PATH');
    oldpath=ENV;
    ENV = strsplit(ENV, ';');
    items_to_add_to_path = {
        char(fullfile(py_root_n2v, 'Library', 'mingw-w64', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'usr', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'bin'))
        char(fullfile(py_root_n2v, 'Scripts'))
        char(py_root_n2v)
        };
    ENV = [items_to_add_to_path(:); ENV(:)];
    ENV = unique(ENV, 'stable');
    ENV = strjoin(ENV, ';');
    
    fprintf('Starting n2v in python ...')
    
    setenv('PATH', ENV);
    [status,cmdout]=system(['python ',fname],'-echo');
    setenv('PATH',oldpath);
    if status==0
        fprintf('N2v finished successfully')
    else
        warning('N2v has a warning. Please check cmdout.')
    end

end

function alignBC2gene(BC_refch,gene_refch,BC_name,gene_name)
% align BC to genes

folders=get_folders();
mkdir('BCtogenealignmentcomparison');
tic
parfor i=1:length(folders)
    cd(folders{i});
    % if genes have been registered, move the original gene files out.
    if isfolder('original')
        cd original
        if ~isempty(dir(['*',gene_name,'*.tif']))
            movefile(['*',gene_name,'*.tif'], '../');
        end
        regfiles=dir('*reg*.tif');
        if ~isempty(regfiles)
            delete *reg*.tif
        end
        cd ..
    end
    %if BC has been registered before, delete the registered files.
    regfiles=dir('*reg*.tif')
    if ~isempty(regfiles)
        delete *reg*.tif
    end
    
    mmalignBCtogene(BC_refch,gene_refch,BC_name,gene_name);
    movefile('comp.tif',['../BCtogenealignmentcomparison/',folders{i},'comp.tif']);
    
    %if genes have been registered, move the original gene files back.
    if isfolder('original')
        movefile(['*',gene_name,'*.tif'], 'original/')
    end
   
    cd ..
end
toc
end


function organize_geneseq()
    % Sort max proj files into sequences.
    fprintf('Putting geneseq files in the correct folders ...')
    genefolders=dir('*geneseq*');
    genefolders(~[genefolders.isdir])=[];
    genefolders=sort_nat({genefolders.name});
    
    mkdir processed
    
    for i=1
        fprintf(genefolders{i})
        cd([genefolders{i}]);
        files=dir('MAX*.tif');
        files=sort_nat({files.name});
        parfor (n=1:length(files),4)
            filename=textscan(files{n},'%s','Delimiter','.');%remove .tif
            filename=filename{1}{1};
            mkdir(['../processed/',filename]);
            copyfile(files{n},['../processed/',filename,'/','geneseq01.tif']);
        end
        cd ..
    end
    
    for i=2:length(genefolders)
        fprintf(genefolders{i})
        cd([genefolders{i}]);
        files=dir('MAX*.tif');
        files=sort_nat({files.name});
        parfor (n=1:length(files),4)
            filename=textscan(files{n},'%s','Delimiter','.');%remove .tif
            filename=filename{1}{1};
            copyfile(files{n},['../processed/',filename,'/','geneseq',num2str(i,'%.2u'),'.tif']);
        end
        cd ..
    end
    fprintf('Done.')
end

function organize_bcseq()
% Sort max proj files into sequences.
fprintf('Putting bcseq files in the correct folders ...')


if ~exist('processed','file')
    mkdir processed;
    isfirstcycle=1;
else
    isfirstcycle=0;
end

%bc seq
bcfolders=dir('*bcseq*');
bcfolders(~[bcfolders.isdir])=[];
bcfolders=sort_nat({bcfolders.name});

if isfirstcycle==1
    for i=1
        fprintf(bcfolders{i})
        cd([bcfolders{i}]);
        files=dir('MAX*.tif');
        files=sort_nat({files.name});
        parfor (n=1:length(files),4)
            filename=textscan(files{n},'%s','Delimiter','.');%remove .tif
            filename=filename{1}{1};
            mkdir(['../processed/',filename]);
            copyfile(files{n},['../processed/',filename,'/','bcseq01.tif']);
        end
        cd ..
    end

    for i=2:length(bcfolders)
        fprintf(bcfolders{i})
        cd([bcfolders{i}]);
        files=dir('MAX*.tif');
        files=sort_nat({files.name});
        parfor (n=1:length(files),4)
            filename=textscan(files{n},'%s','Delimiter','.');%remove .tif
            filename=filename{1}{1};
            copyfile(files{n},['../processed/',filename,'/','bcseq',num2str(i,'%.2u'),'.tif']);
        end
        cd ..
    end
else
    for i=1:length(bcfolders)
        cd([bcfolders{i}]);
        files=dir('MAX*.tif');
        files=sort_nat({files.name});
        for n=1:length(files)
            filename=textscan(files{n},'%s','Delimiter','.');%remove .tif
            filename=filename{1}{1};
            copyfile(files{n},['../processed/',filename,'/','bcseq',num2str(i,'%.2u'),'.tif']);
        end
        cd ..
    end
end
fprintf('Done\n');
end

function cmdout=n2v_geneonly(fname)
    if ~exist('fname','var')
        fname='n2vprocessing.py';
    end

    %This needs to be changed on new systems
    py_root_n2v = fileparts("C:\Users\mara.rue\Anaconda3\envs\n2v_mara\python.exe");
    ENV = getenv('PATH');
    oldpath=ENV;
    ENV = strsplit(ENV, ';');
    items_to_add_to_path = {
        char(fullfile(py_root_n2v, 'Library', 'mingw-w64', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'usr', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'bin'))
        char(fullfile(py_root_n2v, 'Scripts'))
        char(py_root_n2v)
        };
    ENV = [items_to_add_to_path(:); ENV(:)];
    ENV = unique(ENV, 'stable');
    ENV = strjoin(ENV, ';');
    
    fprintf('Starting n2v in python ...')
    
    setenv('PATH', ENV);
    [status,cmdout]=system(['python ',fname],'-echo');
    setenv('PATH',oldpath);
    if status==0
        fprintf('N2v finished successfully')
    else
        warning('N2v has a warning. Please check cmdout.')
    end

end

function cmdout=bardensr()
    %This needs to be changed on new systems
    py_root_n2v = fileparts("C:\Users\mara.rue\Anaconda3\envs\bardensr_mara\python.exe");
    ENV = getenv('PATH');
    oldpath=ENV;
    ENV = strsplit(ENV, ';');
    items_to_add_to_path = {
        char(fullfile(py_root_n2v, 'Library', 'mingw-w64', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'usr', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'bin'))
        char(fullfile(py_root_n2v, 'Scripts'))
        char(py_root_n2v)
        };
    ENV = [items_to_add_to_path(:); ENV(:)];
    ENV = unique(ENV, 'stable');
    ENV = strjoin(ENV, ';');
    
    fprintf('Starting bardensr in python ...')
    cd ..
    setenv('PATH', ENV);
    [status,cmdout]=system('python bardensrbasecall.py','-echo');
    setenv('PATH',oldpath);
    cd processed
    if status==0
        fprintf('bardensr finished successfully')
    else
        warning('bardensr has a warning. Please check cmdout.')
    end

end

function cmdout=cellpose()
    %This needs to be changed on new systems
    py_root_n2v = fileparts("C:\Users\mara.rue\Anaconda3\envs\cell_pose_mara\python.exe");
    ENV = getenv('PATH');
    oldpath=ENV;
    ENV = strsplit(ENV, ';');
    items_to_add_to_path = {
        char(fullfile(py_root_n2v, 'Library', 'mingw-w64', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'usr', 'bin'))
        char(fullfile(py_root_n2v, 'Library', 'bin'))
        char(fullfile(py_root_n2v, 'Scripts'))
        char(py_root_n2v)
        };
    ENV = [items_to_add_to_path(:); ENV(:)];
    ENV = unique(ENV, 'stable');
    ENV = strjoin(ENV, ';');
    
    fprintf('Starting cellpose in python ...')
    cd ..
    setenv('PATH', ENV);
    [status,cmdout]=system('python Cellsegmentation-v065.py','-echo');
    setenv('PATH',oldpath);
    cd processed
    if status==0
        fprintf('Cellpose finished successfully')
    else
        warning('Cellpose has a warning. Please check cmdout.')
    end
end


function register_seq_images(fname,chprofile20x,ball_radius,chshift20x,rgb_intensity)

    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    %worker pools
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",12);
    elseif p.NumWorkers<12
        parpool("local",12);
    end
    fprintf('Starting image registration on %u workers...\n',p.NumWorkers)
    tic
    parfor i=1:length(folders)
    %for i=1
    %    for i=13
        cd(folders{i});

        if isfolder('original')
            cd original
            try
                movefile([fname,'*.tif'],'../')
            catch ME
                %rethrow(ME)
            end
            cd ..
        end
        [~,~,warnmsg]=mmseqalignmentsinglethread_local(fname,chprofile20x,ball_radius,chshift20x);
        if ~isempty(warnmsg)
            warning('%s: %s\n',folders{i},warnmsg);
        end
        cd aligned
        rgbout1(['alignedfixed',fname],rgb_intensity);
        mkdir('../RGB');
        movefile('RGB*.tif','../RGB/');
        cd ../..
    end
    regtime=toc;
    fprintf('Registration finished, total elapsed time is %u hours %u minutes %u seconds',round(regtime/3600),round(rem(regtime,3600)/60),round(rem(regtime,60)))
end

function remake_rgb_img(fname,rgb_intensity)
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    %worker pools
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",24);
    elseif p.NumWorkers<12
        parpool("local",24);
    end
    
    % remake rgb images
    parfor i=1:length(folders)
        %for i=13
        cd(folders{i});
        cd aligned
        rgbout1(['alignedfixed',fname],rgb_intensity);
        mkdir('../RGB');
        movefile('RGB*.tif','../RGB/');
        cd ../..
    end
end


function make_codebook(codebookname)

    load(['../',codebookname],'codebook')
    codebook1=char(codebook(:,2));
    codebookbin=double(size(codebook1));
    codebookbin(codebook1=='G')=8;
    codebookbin(codebook1=='T')=4;
    codebookbin(codebook1=='A')=2;
    codebookbin(codebook1=='C')=1;
    codebookbin=reshape(codebookbin,size(codebook1,1),[]);
    codebookbin=codebookbin*(2.^((4*size(codebook1,2)-4):-4:0))';
    codebookbin=uint8(dec2bin(codebookbin))-48;
    codebookbin=[codebook(:,1),mat2cell(codebookbin,ones(size(codebookbin,1),1))];
    save('codebook.mat','codebook','codebookbin');
    %save codebook for bardensr
    codebookbin1=reshape(cell2mat(codebookbin(:,2)),size(codebookbin,1),4,[]);
    save('../codebookforbardensr.mat','codebookbin1');
end


function organize_hyb_files(hybfolder)

% move files
    if ~exist('hybfolder','var')
        hybfolder='hyb01';
    end
    
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    
    
    cd(['../',hybfolder]);
    
    % copy sequential round images to the correct folders
    for i=1:length(folders)
        copyfile([folders{i},'.tif'],['../processed/',folders{i},'/','nuclear.tif']);
    end
    cd ../processed
end

function register_hyb_images(ch,radius,nuclearch,reg_cycle,chradius,chprofilehyb,chshifthyb, rgbintensity)
    %register hyb to first seq cycle, copy nuclear images
%%  
    clc
    
    p=gcp('nocreate');
    if isempty(p)
        parpool("local",12);
    elseif p.NumWorkers<12
        parpool("local",12);
    end
    
    
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    fprintf(['Registering hyb to seq cycle ', num2str(reg_cycle),' on ', num2str(p.NumWorkers) ' workers...\n'])
    parfor i=1:length(folders)
    %for i=37
        cd(folders{i});
        mmalignhybtoseq('nuclear','geneseq',ch,chprofilehyb,radius,nuclearch,reg_cycle,chradius,chshifthyb);
        rgbout1('aligned',rgbintensity);
        movefile('aligned*.tif','aligned/');
        movefile('RGB*.tif','RGB/');
        cd ..
    
    end
    fprintf('Registration finished, copying files to make a hyb copy ...\n')
    % Make sequential rounds images
    for i=1:length(folders)
        copyfile([folders{i},'/aligned/alignednuclear.tif'],[folders{i},'/aligned/alignedhyb01.tif']);
    end
    fprintf('All done!\n')
end

function pos=stitch_10x_images(ch_count,reg_ch,overlap,x,trim,rescale_factor)
    
    cam_x=3200;
    cam_y=3200;
    fprintf('Stitching whole-slice image.\n');
    niestitchgenessingleplane(ch_count,reg_ch,overlap,x,trim); % This doesn't work reliably with DIC channel, adjust imaging settings.
    
    
    fprintf('Making 10x equivalent stitcjed images.\n')
    % make fake 10x images. 
    mkdir 10x
    cd(fullfile(['stitchchannel',num2str(x)],'ch1stitched'));
    f=dir('*.tif');
    f=sort_nat({f.name});
    for i=1:length(f)
        im=imread(f{i});
        for n=2:5
            im(:,:,n)=imread(['../ch',num2str(n),'stitched/',f{i}]);
        end
        %im=imresize(imrotate(im,180),0.5);
        im=imresize(im,rescale_factor);
        imwrite(im(:,:,1),['../../10x/',f{i}]);
        for n=2:5
            imwrite(im(:,:,n),['../../10x/',f{i}],'WriteMode','Append');
        end
    end
    
    cd ../..
    fprintf('Parse tileconfig for each position.\n');
    % parse tileconfig for each pos, sensitive to file name
%    [folders,pos,xxx,yyy]=get_folders();
            

%     %  Make initial transofrmation
%     itform={};
%     tileconfig={};
%     uniqpos=sort_nat(unique(pos));
%     for m=1:numel(uniqpos)
%         tileconfig{m}=[max(xxx(ismember(pos,uniqpos(m)))),max(yyy(ismember(pos,uniqpos(m))))];
%         itform{m}={};
%         for i=1:tileconfig{m}(2) %y
%             for n=1:tileconfig{m}(1) %x
%                 %itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*2048/2*0.85 (tileconfig{m}(2)-i)*2048/2*0.85 1];%these are theoretical values
%                 itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*(cam_x/2*(1-overlap)) (i-1)*(cam_y/2*(1-overlap)) 1];%these are empirical measures
%             end
%         end
%     end
%     
%     fprintf('Register 20x sequencing images to 10x sequencing images.\n')
%     % Register 40x sequencing images to 10x sequencing images.
%     tic
%     cd 10x
%     allfiles10x=dir('*.tif');
%     allfiles10x=sort_nat({allfiles10x.name});
%     cd ..
%     
%     tform40xto10x={};
%     parfor(i=1:length(folders),4) %not enough memory when slices are large
%         [~,posidx]=ismember(pos{i},uniqpos);
%         file10x=['../../10x/',allfiles10x{posidx}];
%         cd(folders{i})
%         cd original
%         file40x=dir('*seq*.tif');
%         file40x=sort_nat({file40x.name});
%         file40x=file40x{1};
%         itform1=itform{posidx}{yyy(i),xxx(i)};segim1=[];scalediff=rescale_factor;
%         [tform40xto10x{i},~]=mmregister(file40x,file10x,itform1,segim1,scalediff);
%         %copyfile('40xto10xreg.tif',['../../40xto10x/',folders{i},'.tif']);
%         cd ../..
%         
%     end 
%     toc

    %Read 40xto10x tform directly from the ImageJ stitching outputs.
    %reg_ch=5; %for testing
    %x=1;% for testing
    [folders,pos,xxx,yyy]=get_folders();

    fname_all=dir(fullfile(['stitchchannel',num2str(x)],['ch',num2str(reg_ch)],'Pos*.registered.txt'));
    fname_all=sort_nat({fname_all.name})';
    fname_all=cellfun(@(y) fullfile(['stitchchannel',num2str(x)],['ch',num2str(reg_ch)],y),fname_all,'UniformOutput',0);
    
    
    
    %stitchedfname=dir(fullfile('10x','*.tif'));
    %stitchedfname=sort_nat({stitchedfname.name})';
    %stitchedfname=cellfun(@(y) fullfile('10x',y),stitchedfname,'UniformOutput',0);
    %folders=get_folders();
    %
    tforms_converted={};
    
    for i=1:numel(fname_all)
        fname=fname_all{i};
        fid=fopen(fname);
        c=textscan(fid,'%s %*s %s','Delimiter',';');
        fclose(fid);
        tifnames=c{1};
        tforms=c{2};
        I=contains(tifnames,'tif');
        tifnames=tifnames(I);
        tifnames=cellfun(@(x) x(1:end-4),tifnames,'UniformOutput',0);%truncate filenames to get folder names
        tforms=tforms(I);
        tform_xy=cellfun(@(x)textscan(x,'%*f %f %f %*f','Delimiter',',()'),tforms,'UniformOutput',0);
        tform_xy=cell2mat(cellfun(@cell2mat,tform_xy,'UniformOutput',0))/2;
        min_values = min(tform_xy);
        tform_xy = tform_xy+abs(min_values);

    
        
        %stitchedfileinfo=imfinfo(stitchedfname{i});
        %stitched_width=stitchedfileinfo(1).Width;
        tform_xy(:,1)=max(tform_xy(:,1))-tform_xy(:,1);
        tform_xy(:,1)=tform_xy(:,1)-cam_x*trim/2; %add back the trimed portion
        tform_xy(:,2)=tform_xy(:,2)-cam_y*trim/2; %add back the trimed portion
        
        [~,I1]=ismember(tifnames,folders); %match up files
    
    
        for n=1:numel(I1)
            T=[0.5,0,0;0,0.5,0;tform_xy(n,:),1];
            tforms_converted{I1(n)}=affine2d(T);
        end
    
    end
    tform40xto10x=tforms_converted;

    save('40xto10x.mat','tform40xto10x');
    fprintf('All done. Saved transformations to 40xto10x.mat.\n')
end




function check_stitching(tformfname)
    if ~exist('tformfname','var')
        tformfname='40xto10x.mat';
    end
    
    load(tformfname,'tform40xto10x');
    [folders,pos]=get_folders();
    
    mkdir stitched40xto10x
    uniqpos=sort_nat(unique(pos));
    [~,posidx]=ismember(pos,uniqpos);
    %
    parfor i=1:numel(uniqpos)
        %get image size
        %%
        img10xfilename=['10x/',uniqpos{i},'.tif'];
        im10x=double(imread(img10xfilename,1))+double(imread(img10xfilename,2))+double(imread(img10xfilename,3))+double(imread(img10xfilename,4));
        subfolders=find(posidx==i);
        im=zeros(size(im10x));
        %%
        for n=1:numel(subfolders)
            imgfilename=[folders{subfolders(n)},'/original/n2vgeneseq01.tif'];
            im1=double(imread(imgfilename,1))+double(imread(imgfilename,2))+double(imread(imgfilename,3))+double(imread(imgfilename,4));
            im1=imwarp(im1,tform40xto10x{subfolders(n)},'OutputView',imref2d(size(im)));
            im=max(im,im1);
        end
        %%
        imwrite(imresize(uint16(im10x),0.2),['stitched40xto10x/small',uniqpos{i},'.tif']);
        imwrite(imresize(uint16(im),0.2),['stitched40xto10x/small',uniqpos{i},'.tif'],'WriteMode','Append');
        
        
    end
    fprintf('Registered 40x images in /stitched40xto10x. Visually check to make sure they are registered correctly')

end

function [folders,pos,xxx,yyy]=get_folders()
    %get all folder names and corresponding slice names
    folders=dir('MAX*');
    folders=sort_nat({folders.name});
    pos={};
    xxx=[];
    yyy=[];
    for n=1:numel(folders)
        i1=regexp(folders{n},'_');
        pos{n}=folders{n}(i1(1)+1:i1(2)-1);
        xxx(n)=str2double(folders{n}(i1(2)+1:i1(3)-1));
        yyy(n)=str2double(folders{n}(i1(3)+1:end));
    end
    xxx=xxx+1; %xxx and yyy +1
    yyy=yyy+1; %xxx and yyy +1
    
end

function checkregistration(fname,tformfname,scaling)
[folders,pos]=get_folders();
if ~exist('scaling','var')
    scaling=1;
end

if ~exist('tformfname','var')
    tformfname='40xto10x.mat';
end

load(tformfname,'tform40xto10x');
reg_folder_name=['checkregistration'];
mkdir(reg_folder_name)
uniqpos=sort_nat(unique(pos));
[~,posidx]=ismember(pos,uniqpos);
%
f=dir([folders{1},'/RGB/*',fname,'*.tif']);
f={f.name};
parfor i=1:numel(uniqpos)
    %get image size
    %%
    img10xfilename=['10x/',uniqpos{i},'.tif'];
    img10xinfo=imfinfo(img10xfilename);
    subfolders=find(posidx==i);
    im=repmat({zeros(img10xinfo(1).Height,img10xinfo(1).Width,3,'uint8')},numel(f),1);
    %%
    for n=1:numel(subfolders)
        imgfiles=dir([folders{subfolders(n)},'/RGB/*',fname,'*.tif']);
        imgfiles={imgfiles.name};
        for m=1:numel(imgfiles)
            im1=imread([folders{subfolders(n)},'/RGB/',imgfiles{m}]);
            im1=imwarp(im1,tform40xto10x{subfolders(n)},'OutputView',imref2d([size(im{m},1),size(im{m},2)]));
            im1(1:10,1:10,:)=150;%add a border so that it's easy to count
            im{m}=max(im{m},im1);
        end
    end
    
    %%
    for m=1:numel(im)
        im{m}=im{m}*scaling;
        imwrite(imresize(uint8(im{m}),0.2),[reg_folder_name,'/',uniqpos{i},f{m}]);
    end
end
    fprintf('Stitched images in folder: checkregistration%s. Visually check to make sure they are registered correctly.\n',fname)

end


function import_bardensr_results()
    [folders,~]=get_folders();
    load('codebook.mat','codebook');
    id2={};
    lroi={};
    %for i=55
    for i=1:numel(folders)    
        cd(folders{i})
        cd aligned
        try
            m=dlmread('bardensrresult.csv',',',1,0);
            lroi{i}=m(:,[3 2])+1;
            id2{i}=m(:,4)+1;
        catch
            warning('%s has no rolonies\n',folders{i});
            lroi{i}=[];
            id2{i}=[];
        end
    
        
        cd ../..
    end
    save('basecalls.mat','lroi','id2','-v7.3')
    % manually check false positive rate, reiterate basecalling if necessary (keep fpr <5%)
    unusedidx=find(contains(codebook(:,1),'unused'));
    sum(ismember(cell2mat(id2'),unusedidx))
    numel(cell2mat(id2'))
    fpr=sum(ismember(cell2mat(id2'),unusedidx))/numel(cell2mat(id2'))/numel(unusedidx)*sum(~contains(codebook(:,1),'unused'));
    fprintf('Finished importing bardensr results.\n Overall FPR is %.3u\n',fpr);
end

function basecall_hyb(hybthresh,hybbgn)
[folders,~]=get_folders();
    
load('../codebookhyb.mat','codebookhyb');
load('psf.mat','psf')
% read sequential signals.
idhyb={};
lroihyb={};
sighyb={};
parfor(i=1:length(folders),20)
    cd(folders{i});
    cd aligned
    [idhyb{i},lroihyb{i},sighyb{i}]=mmbasecallhyb(codebookhyb,hybthresh,hybbgn,psf);
    cd ../..
end
save('genehyb.mat','idhyb','lroihyb','sighyb','-v7.3');
end

function check_hyb_call(slice)
% manually check the image of a few FOVs
load('../codebookhyb.mat','codebookhyb');
load('genehyb.mat','lroihyb','idhyb');
uniqhybid=unique(idhyb{slice});
figure;
hold on;
for i=1:numel(uniqhybid)
    scatter(lroihyb{slice}(idhyb{slice}==i,1),lroihyb{slice}(idhyb{slice}==i,2),3,'o','filled');
end
pbaspect([1 1 1]);
legend(codebookhyb(:,1),'Location','eastoutside');
set(gca,'ydir','reverse');
end

function import_cellpose(dilation_radius)
if ~exist('dilation_radius','var')
    dilation_radius = 5; %how many pixels to expand on each cell
end
folders=get_folders();
fprintf('Importing Cellpose segmentation results ...\n')

imcell={};%dilated
imcell_original={};%undilated
cellpos={};
celllist={};
parfor i=1:length(folders)
    cd(folders{i});
    cd aligned
    S=load('cellmask.mat');
    maski=S.maski;
    maski_dilated=maski;
    for n=1:dilation_radius
        maski1=imdilate(maski_dilated,ones(3));
        maski_dilated(maski_dilated==0)=maski1(maski_dilated==0);
    end
    %imcell{i}=maski;
    imcell{i}=maski_dilated;
    imcell_original{i}=maski;
    celllist{i}=unique(imcell{i});
    celllist{i}(celllist{i}==0)=[];
    cellpos{i}=zeros(length(celllist{i}),2);
    for n=1:length(celllist{i})
        [y,x]=ind2sub(size(imcell{i}),find(imcell{i}==celllist{i}(n)));
        cellpos{i}(n,:)=[median(x),median(y)];
    end
    cd ../..
end


save('allsegmentation.mat','imcell_original','imcell','cellpos','celllist','-v7.3');
fprintf('Finished saving segmentations to allsegmentation.mat.\n')
end







