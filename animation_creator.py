# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:28:36 2023

@author: Adam Camerer
modified by Barry Kronenfeld

Making videos with ffmpeg (usually comes with windows, maybe other os's also):
 - images must have sequence by numbers (no leading zeroes)
 - first cd into the folder with your image sequence
 - then enter the following command into a command prompt:

ffmpeg -framerate 30 -i k_%d.png -c:v libx264 -r 30 -pix_fmt yuv420p animation.mp4

Parameters you are most likely to want to change:
 10	               input images per second
 k_%d.png          file name; %d will be replaced by numbers
 -r 30             output frames per second
 animation.mp4     output file name
"""
import test_utils as tu
import os
import numpy as np

#set folder and base file path
cases = [(11,0)] # [(0,1), (4,1), (6,0), (8,0), (9,0)]
base_image_folder = r'C:\Users\xsyao\CaGIS Board Dropbox\cantaloupe bob\Barry\Research\Projects\polyline difference metrics\Hausdorff\animation\hausdroff_animation\outputs'

# define animation parameters
dpi = 200 # default image size is 6" x 4"
k_buffer = 0.005 # overlap between successive distance function sections
steps = 300 # number of frames to produce

#wipe outputs folder
if not os.path.isdir(base_image_folder):      
    os.mkdir(base_image_folder)

#loop over all cases
for case_num, seg_num in cases:
    # setup case number folder and file information
    print("case number "+str(case_num))        
    case_num_folder = os.path.join(base_image_folder, "case_"+str(case_num))
    if not os.path.isdir(case_num_folder):              
        os.mkdir(case_num_folder)
    image_file_base = os.path.join(case_num_folder, "k_")

    # get data for given case
    A,B = tu.case(case_num)

    # establish k values for each frame in the animation
    step_size = 1/steps
    # parameterbs_list = tuple(round(k,2) for k in np.arange(0, 1+step_size, step_size))

    def generate_plot(plot_num, steps):
            imagefile = "{}{}.png".format(image_file_base, plot_num)
            
            tu.plot_hausdorff_solution(A, seg_num, B, 
                                show_labels = True, 
                                k=plot_num/steps,
                                k_buffer = k_buffer,
                                imagefile = imagefile,
                                dpi = dpi,
                                k_use_hausdorff=False,
                                stop_distance_function_at_k=True)
    for i in range(steps+1):
        print("generating frame {} of {}".format(i,steps))
        generate_plot(i,steps)
print("finished!!!")