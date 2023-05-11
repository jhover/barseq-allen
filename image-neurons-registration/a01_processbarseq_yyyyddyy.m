%% combine data together
% load('40xto10x_1.mat','tform40xto10x_1');
% load('40xto10x_2.mat','tform40xto10x_2');
% load('40xto10x_3.mat','tform40xto10x_3');
% 
% tform40xto10x=[tform40xto10x_1,tform40xto10x_2,tform40xto10x_3];
% save('40xto10x.mat','tform40xto10x');
% 
% batches=[ones(1,numel(tform40xto10x_1)),ones(1,numel(tform40xto10x_2))*2,3*ones(1,numel(tform40xto10x_3))];
% save('batches.mat','batches');
cd processed
clc
load('40xto10x.mat','tform40xto10x');
batches=ones(1,numel(tform40xto10x)); %Change this based on each batch number! 
save('batches.mat','batches');


% Bardensr basecall and import
% after running bardensr, read in the results.
%%
clc
load('batches.mat','batches');

use_predefined_thresh=0;
cmdout=bardensr();
%% 

import_bardensr_results(batches);
%%
% call hyb rolonies
% % make hybridization codebook
% %codebookhyb=cell(1);%manually copy codebook here. First column is gene name and second column is a mx1 vector for m rounds. only one element is a non-zero value indicating the channel index that correspond to the gene.
% %openvar codebookhyb
% %save('codebookhyb.mat','codebookhyb');
hybthresh=[50 10 40 10];
hybbgn=[50 60 40 30]; %second channel is higher to remove bleedthrough from the txred channel
basecall_hyb(hybthresh,hybbgn);

%
slice=1;
check_hyb_call(slice)

%

%% Segmentation with cellpose

% run cellpose in python, dilate each cell by 3 pixels
cmdout=cellpose();
dilation_radius=3;
import_cellpose(dilation_radius);


%% assign rolonies to cells, transform all coordiates to 10x

assign_gene_rolonies();

rolonies_to_10x();
%
use_mock=1;
calculate_depth(use_mock);


% organize data into common format

startingsliceidx=1; %Change this each batch, SEQUENTIAL for whole brain
data_fname=organize_processed_data(startingsliceidx);
%% 

%
% find overlapping cells and remove them, also set neurons.barcoded if with barcodes.

boxsize=5; % box size for identifying overlaps
pixelsize=6.5/10; % 10x pixel size
filter_overlapping_neurons(data_fname,boxsize,pixelsize);
%% 

rolthresh = [30 30 30 30];
gaussrad=0;
basecall_barcodes(rolthresh, gaussrad); %change how barcodes are "called" - but not for all of them just those in cell


% assign barcodes rolonies to neurons
assign_bc2cell(); %%This is the code I probably have to change

%transform barcodes to 10x coordinates
bc_to_10x();

% determine barcoded neurons by counts

count_thresh=3;
err_corr_thresh=2;% within a cell body, how to collapse barcodes
data_fname='alldata20230105';
mismatch_thresh=1; % allow this mismatch when matching barcodes

organize_bc_data(count_thresh, err_corr_thresh,data_fname,mismatch_thresh);
%% 

% basecall somas and add soma bc data


basecall_somas();

%%
add_somabc();
%%

filter_somabc('filt_neurons-bc-somas.mat', ...
    'complexity',-0.9, ...
    'score',0.85, ...
    'signal',150 ...
    );


%% compress RGB to jpeg, and compress aligned and original to zip
folders=get_folders();
%compress RGB to jpeg
parfor i=1:numel(folders)
    imfiles=dir(fullfile(folders{i},'RGB','*.tif'));
    imfiles={imfiles.name};
    for n=1:numel(imfiles)
        im=imread(fullfile(folders{i},'RGB',imfiles{n}));
        imwrite(im,fullfile(folders{i},'RGB',[imfiles{n}(1:end-3),'jpg']));
    end
    delete(fullfile(folders{i},'RGB','*.tif'));
end
%
folders=get_folders();
parfor i=1:numel(folders)
    if exist(fullfile(folders{i},'original'),'dir')
        zip(fullfile(folders{i},'original.zip'),fullfile(folders{i},'original'));
        rmdir(fullfile(folders{i},'original'),'s');
    else
        fprintf('Folder %s already compressed, skip.\n',fullfile(folders{i},'original'));
    end
    if exist(fullfile(folders{i},'aligned'),'dir')
        zip(fullfile(folders{i},'aligned.zip'),fullfile(folders{i},'aligned'));
        rmdir(fullfile(folders{i},'aligned'),'s');
    else
        fprintf('Folder %s already compressed, skip.\n',fullfile(folders{i},'aligned')); 
    end
end




% cleanup original files
folders=get_folders();
for i=1:numel(folders)
    cd(folders{i});
    delete *.tif
    cd ..
end



%%
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

function cmdout=bardensr(use_predefined_thresh)
    %This needs to be changed on new systems
    if ~exist('use_predefined_thresh','var')
        use_predefined_thresh=0;
    end



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
    if use_predefined_thresh
        [status,cmdout]=system('python bardensrbasecall_predefinedthresh.py','-echo');
    else
        [status,cmdout]=system('python bardensrbasecall.py','-echo');
    end
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
        parpool("local",32);
    elseif p.NumWorkers<12
        parpool("local",32);
    end
    fprintf('Starting image registration on %u workers...\n',p.NumWorkers)
    tic
    parfor i=1:length(folders)
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
        [~,~,warnmsg]=mmseqalignmentsinglethread(fname,chprofile20x,ball_radius,chshift20x);
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
        parpool("local",32);
    elseif p.NumWorkers<12
        parpool("local",32);
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
        parpool("local",32);
    elseif p.NumWorkers<12
        parpool("local",32);
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
    cd stitchchannel1/ch1stitched
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
    [folders,pos,xxx,yyy]=get_folders();
            

    %  Make initial transofrmation
    itform={};
    tileconfig={};
    uniqpos=sort_nat(unique(pos));
    for m=1:numel(uniqpos)
        tileconfig{m}=[max(xxx(ismember(pos,uniqpos(m)))),max(yyy(ismember(pos,uniqpos(m))))];
        itform{m}={};
        for i=1:tileconfig{m}(2) %y
            for n=1:tileconfig{m}(1) %x
                %itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*2048/2*0.85 (tileconfig{m}(2)-i)*2048/2*0.85 1];%these are theoretical values
                itform{m}{i,n}=[1 0 0; 0 1 0; (tileconfig{m}(1)-n)*(cam_x/2*(1-overlap)) (i-1)*(cam_y/2*(1-overlap)) 1];%these are empirical measures
            end
        end
    end
    
    fprintf('Register 20x sequencing images to 10x sequencing images.\n')
    % Register 40x sequencing images to 10x sequencing images.
    tic
    cd 10x
    allfiles10x=dir('*.tif');
    allfiles10x=sort_nat({allfiles10x.name});
    cd ..
    
    tform40xto10x={};
    parfor(i=1:length(folders),4) %not enough memory when slices are large
        [~,posidx]=ismember(pos{i},uniqpos);
        file10x=['../../10x/',allfiles10x{posidx}];
        cd(folders{i})
        cd original
        file40x=dir('*seq*.tif');
        file40x=sort_nat({file40x.name});
        file40x=file40x{1};
        itform1=itform{posidx}{yyy(i),xxx(i)};segim1=[];scalediff=rescale_factor;
        [tform40xto10x{i},~]=mmregister(file40x,file10x,itform1,segim1,scalediff);
        %copyfile('40xto10xreg.tif',['../../40xto10x/',folders{i},'.tif']);
        cd ../..
        
    end 
    toc
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
    fprintf('Stitched images in folder: checkregistration. Visually check to make sure they are registered correctly.\n')

end


function import_bardensr_results(batches)
    [folders,~]=get_folders();
    load('codebookM1all.mat','codebook');
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
    %sum(ismember(cell2mat(id2'),unusedidx))
    %numel(cell2mat(id2'))
    fpr=sum(ismember(cell2mat(id2'),unusedidx))/numel(cell2mat(id2'))/numel(unusedidx)*sum(~contains(codebook(:,1),'unused'));
    fprintf('Finished importing bardensr results.\n Overall FPR is %.3u\n',fpr);
    
    if exist('batches','var')
        uniq_batches=unique(batches);
        num_batches=numel(uniq_batches);
        for i=1:num_batches
            fpr=sum(ismember(cell2mat(id2(batches==uniq_batches(i))'),unusedidx))/numel(cell2mat(id2(batches==uniq_batches(i))'))/numel(unusedidx)*sum(~contains(codebook(:,1),'unused'));
            fprintf('Batch %u, FPR is %.3u\n',uniq_batches(i),fpr);
        end
    end
    
    
            



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

function assign_gene_rolonies()
    load('basecalls.mat','lroi');
    load('allsegmentation.mat','imcell');
    load('genehyb.mat','lroihyb');
    
    folders=get_folders();
    
    cellid={};cellidhyb={};
    parfor i=1:length(folders)
        if ~isempty(lroi{i})
            cellid{i}=mmassignrol2cell(lroi{i},imcell{i});
        else
            cellid{i}=[];
        end
        if ~isempty(lroihyb{i})
            cellidhyb{i}=mmassignrol2cell(lroihyb{i},imcell{i});
        else
            cellidhyb{i}=[];
        end
    end
    %
    save('cellid.mat','cellid','cellidhyb','-v7.3');
end

function rolonies_to_10x()
    load('basecalls.mat','lroi');
    load('genehyb.mat','lroihyb');
    load('allsegmentation.mat','cellpos');
    load('40xto10x.mat','tform40xto10x');
    
    % transform rolony and cell coordinates to 10x
    lroi10x=lroi;
    lroihyb10x=lroihyb;
    cellpos10x=cellpos;
    for i=1:length(lroi)
        if ~isempty(lroi{i})
            lroi10x{i}=transformPointsForward(tform40xto10x{i},lroi{i});
        end
        if ~isempty(lroihyb{i})
            lroihyb10x{i}=transformPointsForward(tform40xto10x{i},lroihyb{i});
        end
        if ~isempty(cellpos{i})
            cellpos10x{i}=transformPointsForward(tform40xto10x{i},cellpos{i});
        end
    end
    save('lroi10x.mat','lroi10x','lroihyb10x','cellpos10x','-v7.3');
end



function calculate_depth(use_mock)

if ~exist('use_mock','var')
    use_mock=0;
end


[folders,pos]=get_folders();
if use_mock==0
    %% mark cortex contour
    %mark the contours of the cortex manually
    cd 10x
    allfiles10x=dir('*.tif');
    allfiles10x=sort_nat({allfiles10x.name});
    contours=mmsetsurface(allfiles10x,4);
    cd ..
    save('contours.mat','contours');
else
    contours=repmat({[]},numel(folders),1);
end

% calculate depth of rolonies and cells (leftover0
% we don't use depths and angle anymore, but keeping them here for
% compatibility with old data format

pixelsize=6.5/10; %10x on the kinetix
depths={};angle={};

load('lroi10x.mat','lroi10x','lroihyb10x','cellpos10x');



%uniqpos=sort_nat(unique(pos));
for i=1:length(folders)
    %[~,posidx]=ismember(pos{i},uniqpos);
    cd(folders{i});
    cd aligned %change folder to save output .mat file from mmcalculatedepthcells to the right place
    [depths{i},angle{i}]=mmcalculatedepthrols(contours{i},lroi10x{i},pixelsize);
    [depthshyb{i},anglehyb{i}]=mmcalculatedepthrols(contours{i},lroihyb10x{i},pixelsize);
    [celldepths{i},cellangle{i}]=mmcalculatedepthrols(contours{i},cellpos10x{i},pixelsize);
    cd ../..
end
save('depths.mat','depths','angle','depthshyb','anglehyb','celldepths','cellangle','-v7.3');
end

function fname=organize_processed_data(startingsliceidx,fname)
if ~exist('fname','var')
    fname=['alldata',char(datetime('today','Format','yyyyMMdd')),'.mat'];
end
[folders,pos]=get_folders();
%
% concatenate lroi10x and id22
%this needs to be more efficient and easier to read.
if ~exist('startingsliceidx','var')
    startingsliceidx=1;%what is the slice number for the first slice. Useful for datasets combining multiple sequencing runs
end
startingfovidx=1;
tilesize=[3200 3200];

idall=[];
idall2=[];
lroi10xall=[];
lroi40xall=[];
depthsall=[];
angleall=[];
lroihyb10xall=[];
lroihyb40xall=[];
idhyball=[];
depthshyball=[];
anglehyball=[];
cellidall=[];
cellidhyball=[];
celllistall=[];
cellpos10xall=[];
cellpos40xall=[];
celldepthsall=[];
cellangleall=[];    
uniqpos=sort_nat(unique(pos));


load('basecalls.mat','id2','lroi');
load('lroi10x.mat','lroi10x','lroihyb10x','cellpos10x');
load('cellid.mat','cellid','cellidhyb');
load('depths.mat','depths','angle','depthshyb','anglehyb','celldepths','cellangle');
load('genehyb.mat','lroihyb','idhyb')
load('allsegmentation.mat','celllist','cellpos')
load('codebookM1all.mat','codebook');
load('../codebookhyb.mat','codebookhyb');
sliceidall=[];
sliceidhyball=[];
cellsliceidall=[];
for i=1:length(folders)
%which slice
    [~,posidx]=ismember(pos{i},uniqpos);
%     sequencing rolonies
    idall2=[idall2;id2{i}];
    lroi10xall=[lroi10xall;lroi10x{i}];
    lroi40xall=[lroi40xall;lroi{i}];
    depthsall=[depthsall;depths{i}];
    angleall=[angleall;angle{i}];
    cellidall=[cellidall;double(cellid{i})+(startingfovidx+i-1)*10000];
    sliceidall=[sliceidall;double(ones(size(lroi{i},1),1))*(posidx+startingsliceidx-1)];
%     sequential rolonies
    idhyball=[idhyball;idhyb{i}];
    lroihyb10xall=[lroihyb10xall;lroihyb10x{i}];
    lroihyb40xall=[lroihyb40xall;lroihyb{i}];
    depthshyball=[depthshyball;depthshyb{i}];
    anglehyball=[anglehyball;anglehyb{i}];
    cellidhyball=[cellidhyball;double(cellidhyb{i})+(startingfovidx+i-1)*10000];
    sliceidhyball=[sliceidhyball;double(ones(size(lroihyb{i},1),1))*(posidx+startingsliceidx-1)];
%   
%     cells
    celllistall=[celllistall;double(celllist{i})+(startingfovidx+i-1)*10000];
    cellpos10xall=[cellpos10xall;cellpos10x{i}];
    cellpos40xall=[cellpos40xall;cellpos{i}];
    
    celldepthsall=[celldepthsall;celldepths{i}];
    cellangleall=[cellangleall;cellangle{i}];
    cellsliceidall=[cellsliceidall;double(ones(numel(celllist{i}),1))*(posidx+startingsliceidx-1)];
    
end

% combine two codebooks and filter out fake rolonies
idhyball1=idhyball+size(codebook,1);
cid=[idall2(:,1);idhyball1];
croi=[lroi10xall;lroihyb10xall];
croi40x=[lroi40xall;lroihyb40xall];
cdepths=[depthsall;depthshyball];
cangles=[angleall;anglehyball];
ccodes=[codebook;codebookhyb];
ccellid=[cellidall;cellidhyball];
csliceidall=[sliceidall;sliceidhyball];
ccelldepths=celldepthsall;
ccellangle=cellangleall;
ccellsliceid=cellsliceidall;

ccelllist=celllistall;
ccellpos=cellpos10xall;

% filter out uncalled rolonies
croi(cid==0,:)=[];
croi40x(cid==0,:)=[];
cdepths(cid==0,:)=[];
cangles(cid==0,:)=[];
ccellid(cid==0,:)=[];
csliceidall(cid==0,:)=[];
cid(cid==0,:)=[];
%

%filter out rolonies that were not called or rolonies at the edge 7% of
%images
edge_filter_size=0.07;%fraction to filter around the edge
I=find(~(cid==0|croi40x(:,1)<tilesize(1)*edge_filter_size|croi40x(:,1)>tilesize(1)*(1-edge_filter_size)| ...
    croi40x(:,2)<tilesize(2)*edge_filter_size|croi40x(:,2)>tilesize(2)*(1-edge_filter_size)));
croi=croi(I,:);
cdepths=cdepths(I,:);
cangles=cangles(I,:);
ccellid=ccellid(I,:);
csliceidall=csliceidall(I,:);
cid=cid(I,:);
croi40x=croi40x(I,:);


% make expression matrix
expmat=sparse(size(ccelllist,1),size(ccodes,1));
cid1=cid(ismember(ccellid,ccelllist));
ccellid1=ccellid(ismember(ccellid,ccelllist));
[~,ccellid1]=ismember(ccellid1,ccelllist);
%
tic
for i=1:size(ccodes,1)

    expmat(:,i)=accumarray(ccellid1,1*(cid1==i),[size(expmat,1),1]);
end
toc




% save data
save processeddata -v7.3

% organize variable names
%fix depths


rolonies=struct;
rolonies.id=cid;
rolonies.pos=croi;
rolonies.pos40x=croi40x;
rolonies.cellid=ccellid;
rolonies.angle=cangles;
rolonies.depth=cdepths;
rolonies.genes=ccodes;
rolonies.slice=csliceidall;

neurons.expmat=expmat;
neurons.id=ccelllist;
neurons.pos=ccellpos;
neurons.depth=ccelldepths;
neurons.angle=ccellangle;
neurons.slice=ccellsliceid;
neurons.pos40x=cellpos40xall;
neurons.genes=ccodes;

save(fname,'rolonies','neurons','-v7.3');
end


function filt_neurons=filter_overlapping_neurons(fname_neurons,boxsize,pixelsize)
load(fname_neurons,'neurons');


if ~exist('boxsize','var')
    boxsize=5;
end

if ~exist('pixelsize','var')
    pixelsize=6.5/10;
end

%
slice=unique(neurons.slice)';

removecells_all=zeros(numel(neurons.id),1);
%
for i=slice
    %
    inslice=neurons.slice==i;
    sliceneuron_fovnum=floor(neurons.id(inslice,:)/10000);
    sliceneuron_pos=neurons.pos(inslice,:);
    sliceneuron_expmat=neurons.expmat(inslice,:);

    % "sort and sweep' to identify overlapping cells
    % make every cell into a 10um square box. This is smaller than a lot of
    % cells, but sincewe are looking for complete overlaps, smaller shoudl
    % be fine.

    sliceneuron_x_l=sliceneuron_pos(:,1)*pixelsize-boxsize/2;
    sliceneuron_x_u=sliceneuron_pos(:,1)*pixelsize+boxsize/2;
    sliceneuron_y_l=sliceneuron_pos(:,2)*pixelsize-boxsize/2;
    sliceneuron_y_u=sliceneuron_pos(:,2)*pixelsize+boxsize/2;
    %sort intervals
    [interval_list_x,AABB_x]=sort([sliceneuron_x_l;sliceneuron_x_u]);
    [~,AABB_x]=sort(AABB_x);
    AABB_x=reshape(AABB_x,[],2);
    [interval_list_y,AABB_y]=sort([sliceneuron_y_l;sliceneuron_y_u]);
    [~,AABB_y]=sort(AABB_y);
    AABB_y=reshape(AABB_y,[],2);
    %for each upper bound, find everything smaller than it in the interval
    %list
    %
    c={};
    c1={};
    %tic
    for n=1:size(AABB_x,1)
        %c{n}=find(AABB_x(n,2)>AABB_x(:,1)&AABB_x(n,2)<AABB_x(:,2)& ...
        %    AABB_y(n,2)>AABB_y(:,1)&AABB_y(n,2)<AABB_y(:,2)& ...
        %    sliceneuron_fovnum~=sliceneuron_fovnum(n)); % cells that are potentially overlapping and from different FOVs
        c1{n}=find(AABB_x(n,2)>AABB_x(:,1)&AABB_x(n,2)<AABB_x(:,2)& ...
            AABB_y(n,2)>AABB_y(:,1)&AABB_y(n,2)<AABB_y(:,2)); %just overlapping
        if ~isempty(c1{n})
            c{n}=c1{n}(sliceneuron_fovnum(c1{n})~=sliceneuron_fovnum(n));
        else
            c{n}=c1{n};
        end
    end
    %numel(cell2mat(c'))
    %numel(cell2mat(c1'))
    fprintf('Slice %u: %u out of %u filtered cells are from different fovs\n',i,numel(cell2mat(c')),numel(cell2mat(c1')));

    %toc



    % Ideally we want to evaluate  cell overlap pixel by pixel, but this may  be unnecesasry.
    % Go through the potentially overlapping list, remove cells with lower read counts.
    %tic

    d=find(~cellfun(@isempty,c));
    remove_neuron=zeros(numel(c),1);
    for m=d
        if remove_neuron(m)==0
            idx=[m;c{m}];
            if sum(idx>numel(remove_neuron))>0
                m
            end
            [~,keep_idx]=max(sum(sliceneuron_expmat(idx,:),2));
            idx(keep_idx)=[];
            remove_neuron(idx)=1;
        end
    end

    %toc
    removecells_all(inslice)=remove_neuron;
 
end

% remove the labeled cells
filt_neurons=neurons;
filt_neurons.expmat=filt_neurons.expmat(removecells_all==0,:);
filt_neurons.id=filt_neurons.id(removecells_all==0,:);
filt_neurons.pos=filt_neurons.pos(removecells_all==0,:);
filt_neurons.depth=filt_neurons.depth(removecells_all==0,:);
filt_neurons.angle=filt_neurons.angle(removecells_all==0,:);
filt_neurons.slice=filt_neurons.slice(removecells_all==0,:);
filt_neurons.pos40x=filt_neurons.pos40x(removecells_all==0,:);

%only for barcoded data
if isfield(filt_neurons,'dom_bc')
    filt_neurons.dom_bc=filt_neurons.dom_bc(removecells_all==0,:);
    filt_neurons.barcoded=cellfun(@(x) ~isempty(x),filt_neurons.dom_bc);

end
if isfield(filt_neurons,'dom_bc_count')
    filt_neurons.dom_bc_count=filt_neurons.dom_bc_count(removecells_all==0,:);
end
if isfield(filt_neurons,'all_bc')
    filt_neurons.all_bc=filt_neurons.all_bc(removecells_all==0,:);
end
if isfield(filt_neurons,'all_bc_count')
    filt_neurons.all_bc_count=filt_neurons.all_bc_count(removecells_all==0,:);
end

save('filt_neurons.mat','filt_neurons','removecells_all','-v7.3');
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
        cd ..bclroi
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

function [bcseq, bcseqC, bcscore, bcsig, bclroi,bcint, bcqual]=basecall_barcodes(rolthresh, gaussrad)
    folders=get_folders();
    parfor i=1:numel(folders)
        cd(folders{i})
        cd aligned
        [bcseq{i},bcseqC{i},bcscore{i},bcsig{i},bclroi{i}]=mmbasecallsinglerol(rolthresh,gaussrad);
        bcint{i}=max(bcsig{i},[],3);
        bcqual{i}=bcint{i}./sqrt(sum(bcsig{i}.^2,3));
        cd ../..
    end
    
    %
    save('bc.mat','bcseq','bcsig','bcqual','bcint','bcseqC','bcscore','bclroi','-v7.3');
    fprintf('Basecalled barcodes.\n')
end


function assign_bc2cell()
    load('bc.mat','bclroi');
    load('allsegmentation.mat','imcell');
    % assign rolonies to cells
    folders=get_folders();
    bccellid={};
    parfor i=1:length(folders)
        if ~isempty(bclroi{i})
            bccellid{i}=mmassignrol2cell(bclroi{i},imcell{i}); %rolonies outside of cells are 0
        else
            bccellid{i}=[];
        end
    end
    %
    save('bccellid.mat','bccellid','-v7.3');
    fprintf('Assigned barcodes to cells.\n')
end

function bc_to_10x()
    % transform rolony coordinates to 10x
    load('bc.mat','bclroi');
    load('40xto10x.mat','tform40xto10x');
    bclroi10x=bclroi;
    for i=1:length(bclroi)
        if ~isempty(bclroi{i})
            bclroi10x{i}=transformPointsForward(tform40xto10x{i},bclroi{i});
        end
    
    end
    save('bclroi10x.mat','bclroi10x','-v7.3');
    fprintf('Transformed barcodes to stitched 10x coordinates.\n')
end


function organize_bc_data(count_thresh,err_corr_thresh,data_fname,mismatch_thresh)
[folders,pos]=get_folders();
dom_bc=cell(numel(folders),1);
all_bc=cell(numel(folders),1);
all_bc_count=cell(numel(folders),1);
dom_bc_count=cell(numel(folders),1);

load('allsegmentation.mat','celllist');
load('bccellid.mat','bccellid');
load('bclroi10x.mat','bclroi10x');
load('bc.mat','bcseq','bcscore','bcint','bcsig','bclroi')


parfor i=1:length(folders)
    
    dom_bc{i}={};
    dom_bc_count{i}=zeros(numel(celllist{i}),1);
    all_bc{i}={};
    all_bc_count{i}={};

    for n=1:numel(celllist{i})

        if sum(bccellid{i}==celllist{i}(n))>count_thresh % if more than threshold number of bc, then count bc molecules
            
            bc=bcseq{i}(bccellid{i}==celllist{i}(n),:); %all barcode sequences in that cell.
            %collapse simialr barcodes.

            [uniqbc,~,ic]=unique(bc,'rows');
            c1=histcounts(ic,[0;unique(ic)]+0.5);
            [sortedc1,I]=sort(c1,'descend');
            uniqbc=uniqbc(I,:); %sorted by abundunce within the cell

            bclist=uniqbc(1,:);bccount=sortedc1(1);
            
            if size(uniqbc,1)>0
                for m=2:size(uniqbc,1)
                    [min_d,I]=min(pdist2(bclist,uniqbc(m,:),'hamming')*size(uniqbc,2));
                    if min_d<=err_corr_thresh
                        bccount(I)=bccount(I)+sortedc1(m);
                    else
                        bclist(end+1,:)=uniqbc(m,:);
                        bccount(end+1)=sortedc1(m);
                    end
                end
            end
            
            [~,I]=max(bccount);
            bclist=int8(bclist);

            if bccount(I)>count_thresh
                dom_bc{i}{n}=bclist(I,:);
                dom_bc_count{i}(n)=bccount(I);
                all_bc{i}{n}=bclist;
                all_bc_count{i}{n}=bccount;
            else
                dom_bc{i}{n}=[];
                dom_bc_count{i}(n)=0;
                all_bc{i}{n}=[];
                all_bc_count{i}{n}=0;
            end
        else
                dom_bc{i}{n}=[];
                dom_bc_count{i}(n)=0;
                all_bc{i}{n}=[];
                all_bc_count{i}{n}=0;
        end

    end
end

% endpos40x
save('cellbc.mat','dom_bc','all_bc');

% combine barcode information with neurons
load(data_fname,'rolonies');
load('filt_neurons.mat','filt_neurons');

filt_neurons.dom_bc=cell(numel(filt_neurons.id),1);
filt_neurons.dom_bc_count=zeros(numel(filt_neurons.id),1);
filt_neurons.all_bc=cell(numel(filt_neurons.id),1);
filt_neurons.all_bc_count=cell(numel(filt_neurons.id),1);

for i=1:numel(dom_bc)
    [~,I]=ismember(double(celllist{i})+i*10000,filt_neurons.id);
    filt_neurons.dom_bc(I(I>0))=dom_bc{i}(I>0);
    filt_neurons.dom_bc_count(I(I>0))=dom_bc_count{i}(I>0);
    filt_neurons.all_bc(I(I>0))=all_bc{i}(I>0);
    filt_neurons.all_bc_count(I(I>0))=all_bc_count{i}(I>0);
end


% bc
bc=struct;
bc.pos40x=cell2mat(bclroi');
bc.pos=cell2mat(bclroi10x');
bc.slice=[];
bc.seq=cell2mat(bcseq');
bc.qual=cell2mat(bcscore');
bc.int=cell2mat(bcint');
bc.sig=cell2mat(bcsig');
bc.matched_cell=[];
%

startingsliceidx=min(filt_neurons.slice);

uniqpos=sort_nat(unique(pos));
[~,posidx]=ismember(pos,uniqpos);
bcslice={};
parfor i=1:numel(bclroi)
    bcslice{i}=double(ones(size(bclroi{i},1),1))*(posidx(i)+startingsliceidx-1);
end
bc.slice=cell2mat(bcslice');



% match bc to cells by calculating projection on the ref barcode vector
% and/or mismatches. qual/int thresh are optional


%filt= min(bc.qual)>qual_thresh & ...
%    min(bc.int)>int_thresh;
% 
%all cellular barcodes
all_cell_bc=cell2mat(filt_neurons.dom_bc(~cellfun(@isempty,filt_neurons.dom_bc)));
all_cell_bc_count=filt_neurons.dom_bc_count(~cellfun(@isempty,filt_neurons.dom_bc));
all_cell_bc_cellid=filt_neurons.id(~cellfun(@isempty,filt_neurons.dom_bc));


% sort by counts, then throw away cells with repeated barcodes
[~,I]=sort(all_cell_bc_count,'descend');
all_cell_bc=all_cell_bc(I,:);
all_cell_bc_cellid=all_cell_bc_cellid(I);
%
[all_cell_bc,ia,~]=unique(all_cell_bc,'rows','first');
all_cell_bc_cellid=all_cell_bc_cellid(ia);
%label cells with repeat barcodes
d1=squareform(pdist(all_cell_bc,'hamming'));
d1=(d1+eye(size(d1)))*size(all_cell_bc,2);
repeat_bc=min(d1,[],2)<=mismatch_thresh;
%

% match bc to cells
d=pdist2(bc.seq,all_cell_bc,'hamming')*size(bc.seq,2);
%
[dmin,I]=min(d,[],2);
idx=(dmin<=mismatch_thresh).*(sum(d==dmin,2)==1).*I;
bc.matched_cell=zeros(size(bc.pos40x,1),1);
bc.matched_cell(idx>0)=all_cell_bc_cellid(idx(idx>0));
%
filt_neurons.repeat_bc=repeat_bc;

save([data_fname(1:end-4),'-bc.mat'],'rolonies','filt_neurons','bc','-v7.3');
save('filt_neurons-bc.mat','filt_neurons');

end

function plot_genes_on_fov(rolonies,gene_list,radius)
if ~exist('radius','var')
    radius=4;
end

imagesize=[3200 3200];

%make images of rolonies on each fov
output_dir=fullfile('..','figs',[gene_list{:}]);
mkdir(output_dir);
%gene_list={'RG(B19)'};
% make images with gene_list

uniq_FOV=unique(floor(rolonies.cellid/10000));
[~,gene_idx]=ismember(gene_list,rolonies.genes(:,1));
%colors=hclrainbow(numel(gene_idx)+1);
%colors=colors(1:numel(gene_idx),:);
colors=colormap('lines');
close
colors=colors(1:7,:);
colors=1-(1-colors)/2;%brighten the colors
if numel(gene_idx)>size(colors,1)
    colors=[colors;distinguishable_colors(numel(gene_idx)-size(colors,1),[colors;0 0 0;1 1 1])];
else
    colors=colors(1:numel(gene_idx),:);
end
folders=get_folders();
%%
for ii=1:numel(uniq_FOV)
    im=zeros([imagesize,numel(gene_idx)]);
    in_slice=floor(rolonies.cellid/10000)==uniq_FOV(ii);
    parfor n=1:numel(gene_idx)
        im(:,:,n)=full( ...
            sparse(rolonies.pos40x(in_slice&rolonies.id==gene_idx(n),2), ...
            rolonies.pos40x(in_slice&rolonies.id==gene_idx(n),1), ...
            1,imagesize(2),imagesize(1)) ...
            );
        im(:,:,n)=imdilate(im(:,:,n),offsetstrel('ball',radius,1));
        im(:,:,n)=im(:,:,n)-min(im(:,:,n),[],'all');
    end
    im=reshape(im,[],size(im,3));
    im_rgb=uint8(reshape(min(im*colors,1),imagesize(1),imagesize(2),[])*255);
    imwrite(im_rgb,fullfile(output_dir,[folders{ii},'.jpg']));
end
    
figure;%imshow(im_rgb);title(folders(ii),'Interpreter','none')
hold on;
for n=1:numel(gene_idx)
    scatter([],[],1,colors(n,:),'filled');
end
legend(gene_list, ...
    'FontSize',10, ...
    'Location','best')
axis off

exportgraphics(gcf,fullfile(output_dir,'legend.pdf'),'ContentType','vector')
close all;
end




