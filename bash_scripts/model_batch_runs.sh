#!/bin/bash

# runs the model wanted. 

if (( $# != 0 ))
then

    echo "Usage $0,   "
	
else

    data_dir_tot='..'/runs_clean_F8_t7/
    #~ data_dir_tot='..'/sims_for_fig2/
    #~ data_dir_tot='..'/sims_for_fig4/
    #~ data_dir_tot='..'/sims_for_fig5/

 
    if [ ! -d $data_dir_tot ]; then
        mkdir $data_dir_tot
    fi        

 
    script_dir='.' # scripts directory
    program_dir=".."/python_scripts/ # scripts directory
    
    
    curr_dir=`pwd`
    
    declare -A par_id=() # associative array, requires bash 4
    count_tot=0
    rm -rf /users/marchi/.cache/matplotlib/.matplotlib_lock-*


    for mf in  '0.7' '0.75' '1' # this example code cycles through the parameters needed to generate Fig3 B,D
    do


    for Pc in  '-1' #  '1e7'  '-1' # negative indicates linear adsorption profile 
    do


        param_dir_tmp="t7_F8real_mf_${mf}" #
        #~ param_dir_tmp="t4_F8real_mf_${mf}" #
        
        param_dir_d=`echo $param_dir_tmp | sed 's/\./d/g' | sed 's/\s\s*/_/g'` # substitutes dots with d, and spaces and tabs with _
        
        param_dir=`echo $param_dir_d | sed 's/d0\+_/d_/g' | sed 's/d\(0*[123456789][123456789]*\)0\+_/d\1_/g'`
        
        
        #~ param_dir="t4_F9"
        #~ param_dir="t7_F10"
        #~ param_dir="t4_F8_60deg"
        #~ param_dir="t4_F8_vsF10"

        
        data_dir_fin=${data_dir_tot}/${param_dir}/
        
        rm -r $data_dir_fin
        
        cd $program_dir
        
        if [ ! -d $data_dir_fin ]; then
            mkdir $data_dir_fin
        fi
        
        
        input=${data_dir_fin}/'params.txt'
        
        echo "$data_dir_fin"
    

        rm $input
        
        echo "# 1 <mot factor>  2 <P_c>   "  >> $input
        echo "$mf "  " $Pc "  >> $input
        
        
        python run_model.py "$data_dir_fin"  &
        
        
        cd $curr_dir
        
        
        cd $curr_dir
        echo " "
        ((count_tot++))

	done
    done 

fi
