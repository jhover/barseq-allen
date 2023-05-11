function [chprofile,tforms,warnmsg]=mmseqalignmentsinglethread_novatraindata(filename,chprofile, bgnradius,shiftmat)
%fix bleedthrough in all tif images in the current folder and align them.
%This version works for the xlight seq-o-matic. This version requires bleedthrough profile 6/10/2019. This version uses a
%single core. This version preserves non-seq frames and allows the use of
%dic image to pre-align. This version adds the option of shifting channels
%to correct for dichroic thickness. shiftmat is n x 2 for n channels. Each
%row is [Tx Ty].

%check if channels need to be shifted
if ~exist('shiftmat','var')
    shiftmat=zeros(4,2);
end
lastwarn('');

%Find all files
files=dir(fullfile([filename,'*.tif']));
files=sort_nat({files.name});
L=size(files);
fixedfilelist=cell(L(2));
mkdir original;

if files

% Fixes bleedthrough for all files, and store output filename list
parfor k=1:L(2)
    currentfile=files{k};
    [~,fixedfilelist{k}]=fixbleed(currentfile,chprofile,bgnradius,shiftmat);
    movefile(files{k},['original/',files{k}]);
end

% Registration
im=imread(fixedfilelist{1},1);
im(:,:,2)=imread(fixedfilelist{1},2);
imwrite(im(:,:,1),['aligned',fixedfilelist{1}]);
imwrite(im(:,:,2),['aligned',fixedfilelist{1}],'WriteMode','Append');

%copyfile(fixedfilelist{1},['aligned',fixedfilelist{1}]);
parfor k=2:L(2)
    currentfile=fixedfilelist{k};
    tforms{k}=alignseq(currentfile,fixedfilelist{1});
end

% Cleanup files
movefile('aligned*.tif','aligned/');
delete('fixed*tif');

warnmsg=lastwarn;
end


% This subfunction fixes bleed through between channels in illumina seq tiff file with Miseq.
function [chprofile,fixedimage]=fixbleed_Miseq(filename,chprofile,radius,shiftmat)
%read file and median filter
info=imfinfo(filename);
im=[];im1=[];
for i=1:4
    im(:,:,i)=medfilt2(imread(filename,i));
end
if length(info)>4
    for i=5:length(info)
        im1(:,:,i-4)=imread(filename,i);
    end
end

%shift chanels
for i=1:size(shiftmat,1)
    im(:,:,i)=imtranslate(im(:,:,i),shiftmat(i,:));
end
 
%fix n2v artifacts in corners
im(im>60000)=0;%camera is 12 bit, so anything close to the 16-bit range is artifact

%bgn subtraction
%radius=20;
im=im-imopen(im,strel('ball', radius, radius));

% Assign bleed through coefficients and bgn values
%im=reshape(uint16(double(reshape(im,[],4))/chprofile),size(im,1),size(im,2),4);%subtract camera baseline and correct for bleeding

fixedimage=strcat('fixed',filename);
outputfile=fixedimage;
imwrite(uint16(im(:,:,1)),outputfile);
for i=4
    imwrite(uint16(im(:,:,i)), outputfile, 'WriteMode','append');
end
if length(info)>4
    for i=5:length(info)
        imwrite(uint16(im1(:,:,i-4)),outputfile,'WriteMode','append');
    end
end
end









function tform=alignseq(imagename,templatename)
%read images and calculate sum of seq channels. 
info=imfinfo(imagename);
im=zeros(info(1).Height,info(1).Width,length(info));
for n=1:length(info)
    %im(:,:,n)=imwarp(imread(imagename,n),tform,'OutputView',Rfixed);
    im(:,:,n)=imread(imagename,n);
end
alignedim=im;

imagesum=double(sum(im(:,:,1:2),3))./double(max(max(sum(im(:,:,1:2),3))));
template1=imread(templatename,1);
%template2=imread(templatename,2);
%template3=imread(templatename,3);
template4=imread(templatename,2);
templatesum=double(template1+template4)./max(max(double(template1+template4)));

%alignment using ECC
% par.transform = 'translation';
% par.levels = 5;
% par.iterations = 100;
% ransacWarp=iat_ecc(imagesum(200:end-200,200:end-200),templatesum(200:end-200,200:end-200),par);
% [M,N]=size(template4);
% for n=1:length(info)
%     [alignedim(:,:,n),~]=iat_inverse_warping(im(:,:,n),ransacWarp,par.transform,1:N, 1:M);
% end
% tform=[];

%align using imregtform
[optimizer,metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/100; %trial 1 conditions
optimizer.GrowthFactor=1.01;
optimizer.Epsilon=optimizer.Epsilon/10;
optimizer.MaximumIterations=optimizer.MaximumIterations*10;

tform = imregtform(imagesum(200:end-200,200:end-200), templatesum(200:end-200,200:end-200), 'translation', optimizer, metric,'PyramidLevels',5);
Rfixed=imref2d(size(templatesum));
for n=1:length(info)
    alignedim(:,:,n)=imwarp(im(:,:,n),tform,'OutputView',Rfixed);
end

%write aligned image to file.
alignedfile=strcat('aligned',imagename);
imwrite(uint16(alignedim(:,:,1)),alignedfile);
for i=2:size(im,3)
    imwrite(uint16(alignedim(:,:,i)),alignedfile, 'WriteMode','append');
end
end







