function [tform40xto10x,segim10x]=mmregister(file40x,file10x,itform1,segim1,scalediff)
%register 40x sequencing images to 10x sequencing images of the same cycle.
%Produce the transformation tform40xto10x, and transform lroi coordinates
%in 40x to those in 10x (roi10x).
%% read files

%if scalediff is not specified, by default to 0.25 (between 40x and 10x
%diff).
if ~exist('scalediff','var')
    scalediff=0.25;
end

im10x=[];
im40x=[];
for i=1:5
    im40x(:,:,i)=imread(file40x,i);
    im10x(:,:,i)=imread(file10x,i);
    
end


%% scale 40x image
    tformscale=affine2d([scalediff 0 0;0 scalediff 0;0 0 1]);

im40xs=imgaussfilt(imwarp(im40x,tformscale),1);

% %% examine initial images
% figure;imshowpair(sum(im40xs,3),sum(im10x,3));
% 
% figure;imshowpair(im40xs(:,:,5),im10x(:,:,5));

itform0=affine2d(itform1);
im10x1=imwarp(im10x,invert(itform0),'OutputView',imref2d(size(im40xs(:,:,1))));


%% initial registration, translation only
%[optimizer,metric] = imregconfig('multimodal');
%optimizer.InitialRadius = optimizer.InitialRadius/5;
%optimizer.MaximumIterations=300;
%movingRegisteredDefault = imregister(sum(im40xs,3),sum(im10x,3),'translation',optimizer,metric);
%tform1 = imregtform(double(sum(im40xs(:,:,1:4),3)),double(sum(im10x1(:,:,1:4),3)),'translation',optimizer,metric);
%tform1 = imregtform(im40xs(:,:,5),im10x(:,:,5),'translation',optimizer,metric,'InitialTransformation',itform0);
tform1 = imregcorr(double(sum(im40xs(:,:,1:4),3)),double(sum(im10x1(:,:,1:4),3)),'translation','Window',1);


%alignment
% [mpoints,fpoints]=cpselect(double(im40xs(:,:,5))./double(max(max(im40xs(:,:,5)))),double(im10x(:,:,5))./double(max(max(im10x(:,:,5)))),'Wait',true);
% tform1=fitgeotrans(mpoints,fpoints,'nonreflectivesimilarity');


% %% refine registration using the initial transformation, allow rotation (can use similarity too, but difference is ~0.4%.
% [optimizer,metric] = imregconfig('multimodal');
% optimizer.InitialRadius = optimizer.InitialRadius/5;
% optimizer.MaximumIterations=1000;        
% tform2= imregtform(double(sum(im40xs(:,:,1:4),3)),double(sum(im10x(:,:,1:4),3)),'rigid',optimizer,metric,'InitialTransformation',tform1);

 %% examine registration
 moved1=imwarp(im40xs,tform1,'OutputView',imref2d(size(im10x1(:,:,1))));
 %moved2=imwarp(im40xs,tform2,'OutputView',imref2d(size(im10x1(:,:,1))));
 imwrite(uint16(sum(moved1(:,:,1:4),3)),'40xto10xreg.tif');
 imwrite(uint16(sum(im10x1(:,:,1:4),3)),'40xto10xreg.tif','WriteMode','Append');
% 
%figure;imshowpair(sum(moved1(:,:,1:4),3),sum(im10x1(:,:,1:4),3));
% 
 %figure;imshowpair(sum(moved2(:,:,1:4),3),sum(im10x(:,:,1:4),3));
%%
tform40xto10x=tform1;
tform40xto10x.T=tformscale.T*tform1.T*itform0.T;
if ~exist('segim1','var')
    segim10x=[];
elseif isempty(segim1)
    segim10x=[];
elseif ~isempty(segim1)
    segim10x=imwarp(segim1,tform40xto10x,'nearest','OutputView',imref2d(size(im10x(:,:,1))));
end

save('segim10x.mat','segim10x','tform40xto10x');
