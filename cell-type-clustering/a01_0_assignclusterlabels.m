%% Assign c1luster labels
clearvars
clc
load filt_bc_neurons.mat
include_subglu=0;
%% 
match_neuron_clustid(filt_neurons,include_subglu);

%% 

% change variable name
clearvars
load neurons
filt_neurons=neurons;
save('filt_neurons-clust-new.mat','filt_neurons','-v7.3');

%%
function neurons=match_neuron_clustid(neurons,include_subglu)
if ~exist('include_subglu','var')
    include_subglu=0;
end
    neuron_clustid=repmat({'qc_filtered'},numel(neurons.id),1);
    neuron_clustname=repmat({'NA'},size(neuron_clustid));
    neuron_subclass=repmat({'NA'},size(neuron_clustid));
        
    % read cluster labels for cells: whole
    c=table2cell(readtable('analysis/whole/cluster.csv','NumHeaderLines',0));
    ano=readtable('analysis/whole/cluster_annotation.csv');
    [~,I]=ismember(cell2mat(c(:,2)),ano.cluster_id);
    % match cluster labels to cells
    neuron_id=cellfun(@(x) str2double(regexp(x,'(\d*)$','match')),c(:,1));
    [~,inslice]=ismember(neuron_id',neurons.id);
    neuron_clustid(inslice)=ano.cluster_name(I);
    neuron_subclass(inslice)=ano.cluster_name(I);

    % read cluster labels for cells: glu
    c=table2cell(readtable('analysis/glu/cluster.csv','NumHeaderLines',0));
    ano=readtable('analysis/glu/cluster_annotation.csv');
    [~,I]=ismember(cell2mat(c(:,2)),ano.cluster_id);
    % match cluster labels to cells
    neuron_id=cellfun(@(x) str2double(regexp(x,'(\d*)$','match')),c(:,1));
    [~,inslice]=ismember(neuron_id',neurons.id);
    neuron_clustid(inslice)=ano.cluster_name(I);
    neuron_subclass(inslice)=ano.cluster_name(I);

%     % read cluster labels for cells: gaba
    c=table2cell(readtable('analysis/gaba/cluster.csv','NumHeaderLines',0));
    ano=readtable('analysis/gaba/cluster_annotation.csv');
    [~,I]=ismember(cell2mat(c(:,2)),ano.cluster_id);
    % match cluster labels to cells
    neuron_id=cellfun(@(x) str2double(regexp(x,'(\d*)$','match')),c(:,1));
    [~,inslice]=ismember(neuron_id',neurons.id);
    neuron_clustid(inslice)=ano.cluster(I);
    neuron_subclass(inslice)=ano.cluster(I);
    

    % cells that didn't make qc are nan, cells that made qc but are not exc neurons are 'NA'

    if include_subglu>0
        c=table2cell(readtable('analysis/subglu/cluster.csv','NumHeaderLines',0));
        ano=readtable('analysis/subglu/cluster_annotation.csv');
        [~,I]=ismember(c(:,2),ano.cluster_id);
        % match cluster labels to cells
        neuron_id=cellfun(@(x) str2double(regexp(x,'(\d*)$','match')),c(:,1));
        [~,inslice]=ismember(neuron_id',neurons.id);
        neuron_clustid(inslice)=ano.cluster_id(I);
        neuron_clustname(inslice)=ano.cluster(I);
        neuron_subclass(inslice)=ano.subclass(I);
    end

    neurons.clustid=neuron_clustid;
    neurons.subclass=neuron_subclass;
    neurons.clustname=neuron_clustname;
    
    save('neurons.mat','neurons','-v7.3')
end

function filt_neurons=filter_overlapping_neurons(neurons,boxsize,pixelsize)
if ~exist('boxsize','var')
    boxsize=5;
end

if ~exist('pixelsize','var')
    pixelsize=6.5/10;
end

%
slice=unique(neurons.slice)';
clustlist=unique(neurons.clustid);

removecells_all=zeros(numel(neurons.id),1);
%
for i=slice
    %
    inslice=ismember(neurons.clustid,clustlist(1:end-1))&neurons.slice==i;
    [~,sliceneuron_clustid]=ismember(neurons.clustid(inslice),clustlist);
    sliceneuron_fovnum=floor(neurons.id(inslice,:)/10000);
    sliceneuron_pos=neurons.pos(inslice,:);
    sliceneuron_expmat=neurons.expmat(inslice,:);


    %     figure('Position',[100 100 1400 1200]);scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
    %     xl=ylim;yl=ylim;
    %     colormap('parula')
    %     pbaspect([range(xl) range(yl) 1]);
    %     hold on;
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
    fprintf('Slice %u: %u out of %u cells are real\n',i,numel(cell2mat(c')),numel(cell2mat(c1')));

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
    %     %% Visual inspection of results
    %     %plot before and after, shouldn't see any difference. also plot removed
    %     %cells, should be around FOVs
    %     figure('Position',[100 100 2400 1200]);
    %     subplot(1,3,1)
    %     scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
    %     xl=ylim;yl=ylim;
    %     colormap('parula')
    %     pbaspect([range(xl) range(yl) 1]);
    %     title('Before')
    %     subplot(1,3,2);
    %
    %     scatter(sliceneuron_pos(~remove_neuron,1),sliceneuron_pos(~remove_neuron,2),10,sliceneuron_clustid(~remove_neuron),'filled');
    %     xl=ylim;yl=ylim;
    %     colormap('parula')
    %     pbaspect([range(xl) range(yl) 1]);
    %     title('After');
    %     hold on;
    %     subplot(1,3,3);
    %
    %     scatter(sliceneuron_pos(remove_neuron>0,1),sliceneuron_pos(remove_neuron>0,2),10,sliceneuron_clustid(remove_neuron>0),'filled');
    %     xl=ylim;yl=ylim;
    %     colormap('parula')
    %     pbaspect([range(xl) range(yl) 1]);
    %     title('Removed cells');
    %     hold on;
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
filt_neurons.clustid=filt_neurons.clustid(removecells_all==0,:);

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

save('filt_neurons.mat','filt_neurons','removecells_all');
end

function filt_neurons=fix_slices(filt_neurons,removecells_all,do_not_fix)
    % manually draw area that include cells to move, then move draw where to
    % move them to.
    if exist('do_not_fix','var')
        if do_not_fix==1
            filt_neurons.fixedpos=filt_neurons.pos;
        else
            do_not_fix=0;
        end

    else
        do_not_fix=0;
    end

    if do_not_fix==0
        slice=unique(filt_neurons.slice)';
        clustlist=unique(filt_neurons.clustid);


        filt_neurons.fixedpos=filt_neurons.pos;
        for i=slice
            inslice=ismember(filt_neurons.clustid,clustlist(1:end-1))&filt_neurons.slice==i;
            [~,sliceneuron_clustid]=ismember(filt_neurons.clustid(inslice),clustlist);
            sliceneuron_pos=filt_neurons.pos(inslice,:);
            sliceneuron_pos_tformed=sliceneuron_pos;

            %%
            while 1
                close
                figure('Position',[100 100 1400 1200]);scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
                xl=ylim;yl=ylim;
                colormap('parula')
                pbaspect([range(xl) range(yl) 1]);
                hold on;
                originalarea=[];
                title(['Slice ',num2str(i),'. Plot original area, or press enter to skip'])
                while 1
                    [x,y]=myginput(1,'topl');
                    if numel(x)==1
                        originalarea=[originalarea;x y];
                        scatter(x,y,10,'r','filled');
                        if numel(originalarea)>2
                            plot(originalarea(end-1:end,1)',originalarea(end-1:end,2)','r','LineWidth',1);

                        end
                    elseif ~isempty(originalarea)
                        plot([originalarea(1,1) originalarea(end,1)],[originalarea(1,2) originalarea(end,2)],'r','LineWidth',1);
                        break
                    else
                        break
                    end
                end
                if isempty(originalarea)
                    break
                end
                title('Plot transformed area')
                tformedarea=[];
                while 1
                    [x,y]=myginput(1,'topl');
                    if numel(x)==1
                        tformedarea=[tformedarea;x y];
                        scatter(x,y,10,'g','filled');
                        if numel(tformedarea)>2
                            plot(tformedarea(end-1:end,1)',tformedarea(end-1:end,2)','g','LineWidth',1);

                        end
                    else
                        plot([tformedarea(1,1) tformedarea(end,1)],[tformedarea(1,2) tformedarea(end,2)],'g','LineWidth',1);
                        break

                    end
                end
                % find cells in original area.

                [in,on] = inpolygon(sliceneuron_pos(:,1),sliceneuron_pos(:,2),originalarea(:,1),originalarea(:,2));
                in=in|on;
                % find transformation using point pairs
                tform=fitgeotrans(tformedarea,originalarea,'pwl');
                % apply transformation to the in points
                [sliceneuron_pos_tformed(in,1),sliceneuron_pos_tformed(in,2)]=transformPointsInverse(tform,sliceneuron_pos(in,1),sliceneuron_pos(in,2));
                %
                close
                figure('Position',[100 100 1400 1200]);scatter(sliceneuron_pos_tformed(:,1),sliceneuron_pos_tformed(:,2),10,sliceneuron_clustid,'filled');
                xl=ylim;yl=ylim;
                colormap('parula')
                pbaspect([range(xl) range(yl) 1]);
                title('Is the result good (Y/N)?');
                [~,~,button]=ginput(1);
                if button==int8('y')
                    break
                end
            end
            filt_neurons.fixedpos(inslice,:)=sliceneuron_pos_tformed;



        end
    end
    save('filt_neurons.mat','filt_neurons','removecells_all')
end

function visualize_deduplication(neurons,filt_neurons,removecells_all,slicenum)

if ~exist('slicenum','var')
    slicenum=2;
end
%slice=unique(neurons.slice)';
clustlist=unique(neurons.clustid);

inslice=ismember(neurons.clustid,clustlist(1:end-1))&neurons.slice==slicenum;
[~,sliceneuron_clustid]=ismember(neurons.clustid(inslice),clustlist);
sliceneuron_pos=neurons.pos(inslice,:);

figure('Position',[100 100 2400 800]);

subplot(1,4,1)
scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
xl=xlim;yl=ylim;
colormap('parula')
pbaspect([range(xl) range(yl) 1]);
title(['Original slice ',num2str(slicenum)])

inslice=ismember(filt_neurons.clustid,clustlist(1:end-1))&filt_neurons.slice==slicenum;
[~,sliceneuron_clustid]=ismember(filt_neurons.clustid(inslice),clustlist);
sliceneuron_pos=filt_neurons.pos(inslice,:);

subplot(1,4,2)
scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
xl=xlim;yl=ylim;
colormap('parula')
pbaspect([range(xl) range(yl) 1]);
title(['De-duplicated slice ',num2str(slicenum)])

inslice=ismember(neurons.clustid,clustlist(1:end-1))&neurons.slice==slicenum;
[~,sliceneuron_clustid]=ismember(neurons.clustid(inslice&removecells_all),clustlist);
sliceneuron_pos=neurons.pos(inslice&removecells_all,:);

subplot(1,4,3)
scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
xl=xlim;yl=ylim;
colormap('parula')
pbaspect([range(xl) range(yl) 1]);
title(['Removed cells slice ',num2str(slicenum)])


% inslice=ismember(filt_neurons.clustid,clustlist(1:end-1))&filt_neurons.slice==i;
% [~,sliceneuron_clustid]=ismember(filt_neurons.clustid(inslice),clustlist);
% sliceneuron_pos=filt_neurons.fixedpos(inslice,:);

% subplot(1,4,4)
% scatter(sliceneuron_pos(:,1),sliceneuron_pos(:,2),10,sliceneuron_clustid,'filled');
% xl=ylim;yl=ylim;
% colormap('parula')
% pbaspect([range(xl) range(yl) 1]);
% title('Fixed cells')
end

