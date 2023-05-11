clc
clearvars
match_CCF('streamLine_registered_filt_neurons.csv','filt_neurons_D076_4L-clust_CCF.mat') % if no CCF registration, use 0 as input to use original coordinates inplace of CCF
%%
function conn_count_normalize=nornalize_subclass_2(p_region,I2,uniq_presyn_subclass,conn_count,subclass_count_mat,subclass_all,presyn_region)
 conn_count_normalize=conn_count;
 unique_I2=unique(I2)
 for i =1:length(unique_I2)
        index=unique_I2(i);
        if index==0
           continue
        end
        subclass=uniq_presyn_subclass(index);
        
        col=find(ismember(subclass_all,subclass))
        subclass_base_count=sum(subclass_count_mat(:,col));
        conn_count_normalize(index)=conn_count(index)/subclass_base_count;
   end    

end


function conn_count_normalize=nornalize_subclass_1(p_region,I2,uniq_presyn_subclass,conn_count,subclass_count_mat,subclass_all,presyn_region)
 conn_count_normalize=conn_count;
 unique_I2=unique(I2)
 for i =1:length(unique_I2)
        index=unique_I2(i);
        if index==0
           continue
        end
        subclass=uniq_presyn_subclass(index);
        row=ismember(presyn_region,p_region)
        col=find(ismember(subclass_all,subclass))
        subclass_base_count=sum(subclass_count_mat(row,col));
        conn_count_normalize(index)=conn_count(index)/subclass_base_count;
   end    

end



function cmap=hclrainbow(num)

hcl=[linspace(0,360,num)',ones(num,1)*50,ones(num,1)*70];
lch=hcl(:,[3 2 1]);
lch=reshape(lch,1,size(lch,1),[]);
cmap=colorspace('lch->rgb',lch);
cmap=squeeze(cmap);
end

function match_CCF(regfname,neuronfname)
%% use a single input 0 if not registering to CCF
if ~exist('regfname','var')
    regfname='neurons_registered.csv';
end
if ~exist('neuronfname','var')
    neuronfname='filt_neurons.mat';
end
load(neuronfname,'filt_neurons')
CCFneuronfname=[neuronfname(1:end-4),'_CCF.mat'];

if regfname==0
    filt_neurons.CCF=[filt_neurons.slice*20,filt_neurons.pos(:,[2 1])];
else
    % read CCF
    ANO=read_ccf();
    
    %get CCF coorinates
    CCF=readtable(regfname);
    neuronid=CCF.id;
    CCF=[CCF.x_CCF,CCF.y_CCF,CCF.z_CCF];
    
    [~,I]=ismember(filt_neurons.id,neuronid);
    filt_neurons.CCF=CCF(I,:);
    filt_neurons.CCF(filt_neurons.CCF(:)<1)=1;
    
    %find CCF ids and matching area names
    %match CCF v3 annotation
    filt_neurons.CCFano=zeros(size(filt_neurons.CCF,1),1);
    for i=1:numel(filt_neurons.CCFano)
        try
            filt_neurons.CCFano(i)=ANO(round(filt_neurons.CCF(i,1)),round(filt_neurons.CCF(i,2)),round(filt_neurons.CCF(i,3)));
        catch
            warning('some cells are outside of the CCF volume.\n');
        end
    end
    [ccf_id,ccf_parent,ccf_acronym]=parse_CCF_ontology('CCF_ontology.json');
    
    filt_neurons.CCFname=cell(size(filt_neurons.CCF(:,1)));
    filt_neurons.CCFparentname=filt_neurons.CCFname;
    %%
    [~,I]=ismember(filt_neurons.CCFano,ccf_id);

    [~,I1]=ismember(ccf_parent(I),ccf_id);
    filt_neurons.CCFname=ccf_acronym(I);
    filt_neurons.CCFparentname=repmat({'NA'},numel(filt_neurons.CCFname),1);
    filt_neurons.CCFparentname(I1>0)=ccf_acronym(I1(I1>0));
    
end
%%
save(CCFneuronfname,'filt_neurons','-v7.3');

end

function ANO=read_ccf(fname,CCFsize)
if ~exist('fname','var')
    fname='annotation.raw';
end
if ~exist('CCFsize','var')
    CCFsize = [528 320 456];
end

fid = fopen(fname, 'r', 'l' );
ANO = fread( fid, prod(CCFsize), 'uint32' );
fclose( fid );
ANO = reshape(ANO,CCFsize);

end

function [ccf_id,ccf_parent,ccf_acronym]=parse_CCF_ontology(fname)
%read CCF v3 ontology
if ~exist('fname','var')
    fname = 'CCF_ontology.json';
end
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
ccf_ontology = jsondecode(str);
%flatten tree into table with name, id, and parent id.
curr_level=ccf_ontology.msg.children;

ccf_areas=curr_level;
n=1;
while 1
    if numel(curr_level)~=0
        % add the first entry in the list to entries
        ccf_areas(n)=curr_level(1);
        % append children to the end of the list
        curr_level=[curr_level;curr_level(1).children];
        %remove the added entry
        curr_level(1)=[];
        n=n+1;
    else
        break
    end
end
ccf_id=[ccf_areas.id];
ccf_parent=[ccf_areas.parent_structure_id];
ccf_acronym={ccf_areas.acronym}';

ccf_id(end+1)=997;
ccf_parent(end+1)=0;
ccf_acronym{end+1}='NA';

ccf_id(end+1)=0;
ccf_parent(end+1)=0;
ccf_acronym{end+1}='NA';
end




function plot_ccf_slices(filt_neurons,plot_others)
% plot neurons per slice
if ~exist('plot_others','var')
    plot_others=1;
end



clc
cellsize=5;
transparency=0.5;

size(unique(filt_neurons.clustid));
[~,clustidx]=ismember(filt_neurons.clustid,unique(filt_neurons.clustid));
[~,excluded_clustidx]=ismember({'non_Exc','NA','NaN'},unique(filt_neurons.clustid));

cmap=hclrainbow(numel(unique(filt_neurons.clustid)));
cmap(ismember(unique(clustidx),excluded_clustidx),:)=0.85;
mkdir('figs_CCF')
mkdir('figs_CCF/Slices_in_CCF')
for i=unique(filt_neurons.slice)'
    %
    figure('Position',[100 100 2400 1800]);
    if plot_others
        I=filt_neurons.slice==i;
    else
        I=filt_neurons.slice==i&~ismember(unique(clustidx),excluded_clustidx);
    end
    scatter(filt_neurons.CCF(I,3),filt_neurons.CCF(I,2),cellsize, ...
        cmap(clustidx(I),:),'filled','MarkerFaceAlpha',transparency,'MarkerEdgeAlpha',transparency)
    set(gca,'ydir','reverse')
    x=range(xlim);y=range(ylim);
    title('Slice ',num2str(i,'%02u'));
    pbaspect([x y 1])
    
    axis off
    exportgraphics(gcf,['figs_CCF/Slices_in_CCF/Slice',num2str(i,'%02u'),'.pdf'],'ContentType','image')
    close;
end
end
