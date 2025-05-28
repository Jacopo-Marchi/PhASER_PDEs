#!/bin/bash

if (( $# != 1 ))
then

	echo "Usage $0, <input_dir> "
	
else

   # py_script_dir=~/Documents/immunitary_systems_viruses_coevo/viral_immune_coevo/python_scripts/contagion_plots # scripts directory
     
    dir_io=$1 # directory with input files
    
    curr_dir=`pwd`
    #~ cd ${dir_io}/plots
    cd ${dir_io}
    #~ cd plots
    j=0
 
    

	#~ cd frames
	#~ ls -v frames/virs_frame*.png | sed 's/^/file /' | sed 's/$/\nduration 0.1/' > movie_input.txt
    #~ for f in ls -v frame_t_*.png
    #~ do
    #~ echo $f
    #~ echo ${f%%[\..]*}
    #~ mv $f ${f%%[\..]*}.png
    #~ done
    
    #~ newdir="pngs"
    #~ mkdir $newdir
    
    #~ for file in frame_t_*.pdf
    #~ do
      #~ convert -density 150 "$file" -quality 90 "${newdir}/${file%.pdf}.png"
    #~ done
    
    
    #~ cd $newdir    
    

	#~ gls -v frame_t_*.png | sed 's/^/file /' | sed 's/$/\nduration 0.5/' > movie_input.txt
    
    
    gls -v frame_time_*.png | sed 's/^/file /' | sed 's/$/\nduration 0.5/' > movie_input.txt

    
	ffmpeg -f concat -i movie_input.txt -y -c:v libx264 -vf scale=1280:-2 -pix_fmt yuv420p movie_frames.mp4
		
		
	
	#~ cd ${dir_io}/plots
    
		
    cd $curr_dir
   
fi  	
        
    
