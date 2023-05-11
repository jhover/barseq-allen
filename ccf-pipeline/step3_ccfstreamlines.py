#!/usr/bin/env python
# coding: utf-8

# In[6]:


import matplotlib.pyplot as plt
import numpy as np
import json
import ccf_streamlines.projection_removed_singlets as ccfproj
import os
import pandas as pd
import ccf_streamlines.morphology as ccfmorph
import seaborn as sns
import time
with open("avg_layer_depths.json", "r") as f:
    layer_tops = json.load(f)

layer_thicknesses = {
        'Isocortex layer 1': layer_tops['2/3'],
        'Isocortex layer 2/3': layer_tops['4'] - layer_tops['2/3'],
        'Isocortex layer 4': layer_tops['5'] - layer_tops['4'],
        'Isocortex layer 5': layer_tops['6a'] - layer_tops['5'],
        'Isocortex layer 6a': layer_tops['6b'] - layer_tops['6a'],
        'Isocortex layer 6b': layer_tops['wm'] - layer_tops['6b'],
}
bf_boundary_finder = ccfproj.BoundaryFinder(
    projected_atlas_file="flatmap_butterfly.nrrd",
    labels_file="labelDescription_ITKSNAPColor.txt",
)

#We get the left hemisphere region boundaries with the default arguments
bf_left_boundaries = bf_boundary_finder.region_boundaries()

# And we can get the right hemisphere boundaries that match up with
# our projection if we specify the same configuration
bf_right_boundaries = bf_boundary_finder.region_boundaries(
    # we want the right hemisphere boundaries, but located in the right place
    # to plot both hemispheres at the same time
    hemisphere='right_for_both',

    # we also want the hemispheres to be adjacent
    view_space_for_other_hemisphere='flatmap_butterfly',
)
ccf_coord_proj = ccfproj.IsocortexCoordinateProjector(
    projection_file="flatmap_butterfly.h5",
    surface_paths_file="surface_paths_10_v3.h5",
    closest_surface_voxel_reference_file="closest_surface_voxel_lookup.h5",
    layer_thicknesses=layer_thicknesses,
    streamline_layer_thickness_file="cortical_layers_10_v2.h5",
)

neuron_df=pd.read_csv('registered_filt_neurons.csv')
neuron_df_filtered =neuron_df[(neuron_df['x_CCF'] <=528) & (neuron_df['y_CCF'] <=320) & (neuron_df['z_CCF'] <=456) & (neuron_df['x_CCF']>=0)& (neuron_df['y_CCF']>=0)& (neuron_df['z_CCF']>=0) ].copy()
neuron_df_filtered['x_ccf_rescale']=neuron_df_filtered['x_CCF']*25
neuron_df_filtered['y_ccf_rescale']=neuron_df_filtered['y_CCF']*25
neuron_df_filtered['z_ccf_rescale']=neuron_df_filtered['z_CCF']*25

n = 100  #chunk row size
list_df = [neuron_df_filtered[i:i+n] for i in range(0,neuron_df_filtered.shape[0],n)]
list_df = np.array_split(neuron_df_filtered, n)
col=neuron_df_filtered.columns.tolist()+['streamline_dim1','streamline_dim2','streamline_dim3','streamline']
neuron_filtered_streamline_df=pd.DataFrame(columns=col)
for i in range(len(list_df)):
    df=list_df[i]
    coor=df[['x_ccf_rescale','y_ccf_rescale','z_ccf_rescale']].values
    try:
        all_coords_slab = ccf_coord_proj.project_coordinates(
        coor,
        thickness_type='normalized_full',
        drop_voxels_outside_view_streamlines=True,
        view_space_for_other_hemisphere='flatmap_butterfly',)
        df['streamline_dim1']=[i[0] for i in all_coords_slab]
        df['streamline_dim2']=[i[1] for i in all_coords_slab]
        df['streamline_dim3']=[i[2] for i in all_coords_slab]
        df['streamline']='good'
    except:
        df['streamline_dim1']=np.nan
        df['streamline_dim2']=np.nan
        df['streamline_dim3']=np.nan
        df['streamline']='bad'
   
    frames =[neuron_filtered_streamline_df,df]
    neuron_filtered_streamline_df=pd.concat(frames)

neuron_filtered_streamline_df_1=neuron_filtered_streamline_df[['id','streamline_dim1','streamline_dim2','streamline_dim3']].copy()
neuron_df_all=neuron_df.merge(neuron_filtered_streamline_df_1,how='left',on='id')
neuron_df_all.to_csv('streamline_registered_filt_neurons.csv',index=False)




# # check to make sure everything is good
# plot_df=neuron_filtered_streamline_df[~(neuron_filtered_streamline_df['streamline_dim1'].isna())]
# sns.scatterplot(data=plot_df,x='streamline_dim1',y='streamline_dim2',s=0.2)
# for k, boundary_coords in bf_left_boundaries.items():
#     plt.plot(*boundary_coords.T, c="black", lw=0.5)
# for k, boundary_coords in bf_right_boundaries.items():
#     plt.plot(*boundary_coords.T, c="black", lw=0.5)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




