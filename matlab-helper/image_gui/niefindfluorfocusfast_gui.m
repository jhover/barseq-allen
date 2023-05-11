function offset=niefindfluorfocusfast_gui(posfname,fpath,focusfolder,initz,endz,stepz,ch,dicch)
% given a list of positions and a z-stack of fluorescent images, find the
% focus offset for each position.
% this version works on 5x shrinked image to speed up.
%% new version to allow defining which channel to use for sequencing, and to split out dic channel
if ~exist('ch','var')
    ch=1;
end 

if ~exist('dicch','var')
    ch=0;
end  

%%
%%tic
data=readmatrix(posfname);
data2=data;
offset=zeros(size(data,1),1);
%%toc
%%
cd(fpath);
%%
%%tic
files=dir('*.tif');
files=sort_nat({files.name});

fileinfo=imfinfo(files{1});
pzc=cell(length(files),3);
for m=1:length(files)
    pzc(m,:)=textscan(files{m},'%*u %u %u %u %*s','Delimiter',{'xy','z','c','.'});
end
pzc=cell2mat(pzc);
%%toc
%%
scaling = 0.5;
cropsize = 500;
best_z=[];
parfor n=1:size(data,1)
    %%
    %tic
    fidx=find(pzc(:,1)==n&pzc(:,3)==ch);
    im=[];
    for i=1:numel(fidx)
        im1=imread(files{fidx(i)});
        im(:,:,i)=imresize(im1(cropsize:end-cropsize,cropsize:end-cropsize),scaling,'method','bicubic');
    end
    %m=zeros(size(im,3),1);
    %toc
    %%
    %tic
    %calculate focus, using mean intensity excluding the top 0.05%
    %%
    %tic
    im1=reshape(im,[],size(im,3));
    im1=sort(im1,'ascend');
    %im1=sort(im1(1:100:size(im1,1),:),'ascend');%downsize by 100x to speedup
    %im1(im1>prctile(im1,99.95))=[];
    [~,I]=max(mean(im1(1:end-round(size(im1,1)*0.0005),:)));%remove the brightest 0.05% in case there are super bright objects dominating
    %toc
    
    if endz>=initz %the piezo and manual focus are reversed
        offset(n)=initz+(I-1)*stepz;
    else
        offset(n)=initz-(I-1)*stepz;
    end
    best_z(n)=I;
    %toc
end
cd ..
%%
data2(:,3)=data2(:,3)+offset;%assuming same direction betwen the piezo and the focus
writematrix(data2,['offset',posfname],'Delimiter',';')
%% If dic channel is present, split off the in-focus dic image.
if dicch>0
    mkdir(['dic',focusfolder]);
    for n=1:size(data,1)
        fidx=find(pzc(:,1)==n&pzc(:,3)==dicch&pzc(:,2)==best_z(n));
        copyfile([focusfolder,'/',files{fidx}],['dic',focusfolder,'/xy',num2str(n,'%.2u'),'z',num2str(best_z(n),'%.2u'),'c',num2str(dicch,'%.2u'),'.tif']);
    end
end





end