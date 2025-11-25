#!/bin/bash


############################################################
# Script Name: model_batch_runs.sh
# Description: cycles through some parameters and runs the PhASERs model
# Author:      Jacopo Marchi
# Date:        11/18/2025
# Usage:       ./model_batch_runs.sh 
# Notes:       This script is a wrapper calling the python script run_model.py. As an example, the default uncommented valued run the data to generate Fig.2A for phage T7. To produce the data for phage T4, uncomment the corresponding lines. Then Figure 2 can be plotted using the corresponding python plot script (fig2.py in this repository)
############################################################


if (( $# != 0 ))
then

    echo "Usage $0,   "
	
else

    # choose output folder
    #~ data_dir_tot='..'/runs_clean_F8_t7/
    #~ data_dir_tot='..'/runs_F8_t7_new_ingredients/
    data_dir_tot='..'/sims_for_fig2/
    #~ data_dir_tot='..'/sims_for_fig4/
    #~ data_dir_tot='..'/sims_for_fig5/

 
    if [ ! -d $data_dir_tot ]; then
        mkdir $data_dir_tot
    fi        

 
    script_dir='.' # scripts directory
    program_dir=".."/python_scripts/ # scripts directory
    
    
    curr_dir=`pwd`
    
    count_tot=0
    rm -rf /users/marchi/.cache/matplotlib/.matplotlib_lock-*
    
    # CYCLE OVER DESIRED INPUT PARAMETERS. THE EXAMPLE UNCOMMENTED VALUES RUN FIGURE 2A


    #~ for mf in  '0.7' '0.75' '1' # this example code cycles through the parameters needed to generate Fig3 B,D
    for mf in   '1' # this example code only runs the null model where infected cells are as fast as non infected ones
    do

    for Pc in  '-1' #  '1e7'  '-1' # negative indicates linear adsorption profile, for instance for T7
    #~ for Pc in  '1e7'  #  '1e7'  '-1' # use this for T4
    do


    #~ for gr_fact in  '0.1' '0.25' '0.5' '2'  '4'  #   # multiplicative factor to modulate maximum growth rate from the strain-specific inferred one. It modifies the yield
    #~ for gr_fact in '0.25' # '0.25' '0.5' '2'  '4'  #   # multiplicative factor to modulate maximum growth rate from the strain-specific inferred one. It modifies the yield
    for gr_fact in  '1' #   # default
    do

    #~ for phi in  '2e-9' #   # cycle through desired adsorption rates
    for phi in  '2e-8' #   # default T7
    #~ for phi in  '6e-08' #   # default T4
    do
    
    #~ for eta in  '1' '15' #   #  cycle through desired lysis rates
    for eta in  '5' #    # default T7 
    #~ for eta in  '2' #   # default T4
    do

        # name output dir according to chosen parameters
        #~ param_dir_tmp="t7_F8_mf_${mf}" #
        #~ param_dir_tmp="t4_F8_mf_${mf}" #
        
        #~ param_dir_tmp="t4_F8_mf_${mf}_L_1" # L=2 directory, L is set in the python code directly. If modulating L, uncomment and rename according to chosen parameters, and comment default dir name below

        #~ param_dir_tmp="t4_F9_mf_${mf}_gr_fact_${gr_fact}_phi_${phi}_eta_${eta}" # GENERAL DIRECTORY NAME. Use when cycling. Name phage according to saturation Pc above. T7 is negative
        #~ param_dir_tmp="t7_F9_mf_${mf}_gr_fact_${gr_fact}_phi_${phi}_eta_${eta}" #


        #~ param_dir_tmp="t7_F9_mf_${mf}_gr_fact_${gr_fact}_rescyield" # modulating the yield (change gr_fact above)
        

        
        param_dir_d=`echo $param_dir_tmp | sed 's/\./d/g' | sed 's/[[:space:]]\{1,\}/_/g'` # substitutes dots with d, and spaces and tabs with _
        
        param_dir=`echo $param_dir_d | sed 's/d0\+_/d_/g' | sed 's/d\(0*[123456789][123456789]*\)0\+_/d\1_/g'`
        
        
        #~ param_dir="t4_F9" # use this for T4, Fig2B
        param_dir="t7_F9"  # use this for T7, Fig2A
        
        #~ param_dir="t4_F8_60deg" # PhASERs colliders
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
        
        echo "# 1 <mot factor>  2 <P_c>  3 <gr_fac>  4 <phi>  5 <eta>   "  >> $input
        echo "$mf "  " $Pc " " $gr_fact " " $phi " " $eta "  >> $input
        
        
        python run_model.py "$data_dir_fin"  &
        
        
        cd $curr_dir
        
        
        cd $curr_dir
        echo " "
        ((count_tot++))

	done
    done 
    done 
    done 
    done

fi
