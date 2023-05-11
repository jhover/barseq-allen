%% plot color coded neurons in CCF

%%
clearvars
clc
load filt_neurons-clust_D077_2L-allcells_CCF.mat
%%
%reg_dir=fullfile('image_registration','newsliceimages','small');
%output_dir=fullfile('figs','CCF_images');
T=readtable('streamline_registered_filt_neurons.csv');
%%
%CCF=[T.x_CCF,T.y_CCF,T.z_CCF];
streamlines=[T.streamline_dim1,T.streamline_dim2,T.streamline_dim3];
[~,I]=ismember(filt_neurons.id,T.id);
%filt_neurons.CCF=CCF(I,:);
filt_neurons.CCF_streamlines=streamlines(I,:);
save('filt_neurons-bc-clust-CCF.mat','filt_neurons','-v7.3');

%%
load filt_neurons-bc-clust-CCF.mat
%%
size(unique(filt_neurons.clustid));
[~,clustidx]=ismember(filt_neurons.clustid,unique(filt_neurons.clustid));
cmap=hclrainbow(numel(unique(filt_neurons.clustid)));

% 
% figure;
% hold on;
% %for i=numel(filt_neurons.slice)
%     %I=filt_neurons.slice==i&clustidx<118;
%     I=clustidx<118;
%     scatter3(filt_neurons.CCF(I,3),filt_neurons.CCF(I,2),filt_neurons.CCF(I,1)+rand(sum(I),1)*8-4,1, ...
%         cmap(clustidx(I),:),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
% %end
% set(gca,'ydir','reverse')
% pbaspect([1 1 1])
%%
cmap=hclrainbow(numel(unique(filt_neurons.clustid)));
cmap(118,:)=0.85;
v=VideoWriter('vid3dww.avi');
v.Quality=95;
open(v);
pause(1)
figure('Position',[300 0 1600 1400]);
hold on;
for i=1:40
    I=filt_neurons.slice==i&clustidx<119;
    scatter3(-filt_neurons.CCF(I,1)+rand(sum(I),1)*8-4,filt_neurons.CCF(I,3),-filt_neurons.CCF(I,2),10, ...
        cmap(clustidx(I),:),'filled','MarkerEdgeAlpha',0.02,'MarkerfaceAlpha',0.2);
    %set(gca,'ydir','reverse');
    pbaspect([528 456 320])
    view(60,30);
    set(gca,'ylim',[0 456],'zlim',[-320 0],'xlim',[-528,0]);
    set(gcf,'color','w');
    axis off
    axis vis3d
    frame = getframe(gcf);
%     for n=1:10
%         writeVideo(v,frame);
%     end
end
%% 

for n=1:360
    view(60+n,30);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)


%%
% 
% 
% 
% 
% function cmap=hclrainbow(num)
% 
% hcl=[linspace(0,360,num)',ones(num,1)*50,ones(num,1)*70];
% lch=hcl(:,[3 2 1]);
% lch=reshape(lch,1,size(lch,1),[]);
% cmap=colorspace('lch->rgb',lch);
% cmap=squeeze(cmap);
% end


