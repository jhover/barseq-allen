function niestitchgenessingleplane(ch,refch,overlap,filenum,cropfrac,fname)

%stitch 1st cycle sequencing images procesed by by reading positions from filename. Output stitched maxprojection images.
%optional: filenum, which cycle to stitch.

%make max projections of stacks given channel numbers, stitch by reading positions from filename. Output stitched maxprojection images.
%this version does not make the final multi-channel tiff to accomodate
%larger files
%cropfrac: crop this fraction of the image before stitching

%% initialize

%ch=3; %number of channels
%overlap=0.15; %fraction of image overlap
%refch=2; %channel used for stitching
if ~exist('cropfrac','var')
    cropfrac=0;
end
%%
if ~exist('filenum','var')
    filenum=1;
end
if ~exist('fname','var')
    fname='geneseq';
end

%% Make max projections and flip both x and y, also collects file names.
files=dir('MAX*');
files=sort_nat({files.name});
sipxy=cell(length(files),3);




for i=1:length(files)
    sipxy(i,:)=textscan(files{i},'%*s %s %u %u',1,'Delimiter','_');
    
end

xy=cell2mat(sipxy(:,2:3));
pos=cellfun(@cell2mat,sipxy(:,1),'UniformOutput',false);
uniqpos=sort_nat(unique(pos));

%build common filename
c=textscan(files{i},'%s','Delimiter','_');
c=c{1}';
c(2,:)={'_'};
c=cell2mat(c(6:end-1));

if contains(fname,'geneseq')
    foldername=['stitchchannel',num2str(filenum)];
else
    foldername=['stitchchannel',fname,num2str(filenum)];
end

mkdir(foldername);

%% write channel to be used for stitching in stitch folder




for n=1:ch
    mkdir([foldername,'/ch',num2str(n)]);
    for i=1:length(files)
        cd(files{i});
        cd original;
        f=dir(['*',fname,'*.tif']);
        f=sort_nat({f.name});
        %im=imrotate(imread(f{filenum},n),180); %read and invert image.
        im=flip(imread(f{filenum},n),2);%read and flip x axis
        if cropfrac>0
            im=im(floor(1+size(im,1)*cropfrac):floor(size(im,1)*(1-cropfrac)),floor(1+size(im,2)*cropfrac):floor(size(im,2)*(1-cropfrac)));%crop image
        end
        imwrite(im,['../../',foldername,'/ch',num2str(n),'/',files{i},'.tif']);
        cd ../../
    end
end

%%	Stitch reference channel using ImageJ
addpath('C:\Fiji.app\scripts');
%javaaddpath 'C:\Program Files\MATLAB\R2018b\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2022a\java\jar\mij.jar'
%javaaddpath 'C:\Program Files\MATLA\R2018b\java'
Miji(false);
cd([foldername]);

mkdir stitchedref
for i=1:length(uniqpos)
    %for each position, count the number of columns and rols
    colnum=max(xy(ismember(pos,uniqpos(i)),1))+1;
    rownum=max(xy(ismember(pos,uniqpos(i)),2))+1;
    %first stitch using ImageJ
    MIJ.run('Grid/Collection stitching', ['type=[Filename defined position] ',...
        'order=[Defined by filename         ] ',...
        'grid_size_x=',num2str(colnum),' grid_size_y=',num2str(rownum),' '....
        'tile_overlap=',num2str(floor(overlap*100)),' ',...
        'first_file_index_x=0 first_file_index_y=0 ',...
        'directory=[',[pwd,'\ch',num2str(refch)],'] ',...
        'file_names=',['MAX_',uniqpos{i},'_{xxx}_{yyy}','.tif '],...
        'output_textfile_name=',[uniqpos{i},'TileConfiguration.txt '],...
        'fusion_method=[Linear Blending] ',...
        'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ',...
        'compute_overlap ',...
        'computation_parameters=[Save computation time (but use more RAM)] ',...
        'image_output=[Fuse and display]']);
    MIJ.run("Flip Horizontally");
    MIJ.run('Save',['Tiff..., path=[',pwd,'\stitchedref\','stitchedref_',uniqpos{i},'.tif]']);
    MIJ.run("Close All")
end


%% transform other channels.
for n=1:ch
    mkdir(['ch',num2str(n),'stitched']);
    for i=1:length(uniqpos)
        MIJ.run('Grid/Collection stitching', ['type=[Positions from file] ',...
            'order=[Defined by TileConfiguration] ',...
            'directory=[',[pwd,'\ch',num2str(n)],'] ',...
            'layout_file=',['..\ch',num2str(refch),'\',uniqpos{i},'TileConfiguration.registered.txt '],...
            'fusion_method=[Linear Blending] ',...
            'regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ',...
            'computation_parameters=[Save computation time (but use more RAM)] ',...
            'image_output=[Fuse and display]']);
        MIJ.run("Flip Horizontally");
        MIJ.run('Save',['Tiff..., path=[',pwd,'\ch',num2str(n),'stitched\',uniqpos{i},'.tif]']);
        MIJ.run("Close All")
    end
end
MIJ.exit;
% %% piece channels together.
% for i=1:length(uniqpos)
%     im=imread(['ch1stitched/',uniqpos{i},'.tif']);
%     imwrite(im,['stitched',uniqpos{i},'.tif']);
%     for n=2:ch
%         im=imread(['ch',num2str(n),'stitched/',uniqpos{i},'.tif']);
%         imwrite(im,['stitched',uniqpos{i},'.tif'],'WriteMode','Append');
%     end
% end
% 
cd ..



%%







