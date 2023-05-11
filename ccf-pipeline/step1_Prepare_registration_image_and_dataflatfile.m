%% load data
clearvars
clc
load filt_neurons-bc-clust

scaling_factor=0.2;
%% For each slice, pouplate sparse matrix, then imwrite
%convert clustid to numbers
uniqclust=unique(filt_neurons.clustid);
[~,clustidx]=ismember(filt_neurons.clustid,uniqclust);
clustidx=max(clustidx)-clustidx+1;%invert the sequence so that the non-exc and qc filtered gets gray

%convert to colormap
cmap=parula(numel(unique(clustidx)));
cmap=[0.5 0.5 0.5;0.3 0.3 0.3;cmap]; %add gray to the first one

%max x and y coordinates
maxx=ceil(max(filt_neurons.pos(:,1)));
maxy=ceil(max(filt_neurons.pos(:,2)));


%% generate cell images for CCF registration
output_dir=fullfile('newsliceimages','small');
mkdir(output_dir);

%mkdir newsliceimages
%cd newsliceimages
%mkdir small
for i=unique(filt_neurons.slice)'
    
   %make sparse matrix
    in_slice=filt_neurons.slice==i;%indices of neurons within this slice
    pos_coor=filt_neurons.pos(:,1)>1&filt_neurons.pos(:,2)>1;
    img={};
    for n=1:3
        %generate sparse matrix
        img{n}=sparse(round(filt_neurons.pos(in_slice&pos_coor,2)), ...
            round(filt_neurons.pos(in_slice&pos_coor,1)), ...
            cmap(clustidx(in_slice&pos_coor),n), ...
            maxy, ...
            maxx ...
            );
        %convolve to make dots bigger
        img{n}=conv2(full(img{n}),ones(20),'same');
        % cap at max 1
        img{n}=min(img{n},1);
    end
    
    %imwrite(cat(3,img{1},img{2},img{3}),['Slice',num2str(i,'%.2u'),'.jpg']);
    imwrite(imresize(cat(3,img{1},img{2},img{3}),scaling_factor),fullfile(output_dir,['Slice',num2str(i,'%.2u'),'.jpg']));
    
    
    
end
%% write filt_neuron coordinates in csv
T=table(filt_neurons.id, ...
    filt_neurons.pos(:,1), ...
    filt_neurons.pos(:,2), ...
    filt_neurons.slice, ...
    filt_neurons.subclass, ...
    filt_neurons.clustid, ...
    'VariableNames', ...
    {'id','x','y','slice','subclass','clustid'} ...
    );



filename='filt_neurons.csv';
writetable(T,filename)