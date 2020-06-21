bool master(char in_traj_PATH[], char out_traj_PATH[], int filecount){

    char *filler = "-----------------------";

    bool   should_open_file;
    bool caught_nan = false;

    int    stat;
    int    count = 0, iter = 0, inrange_iter = 0, file_number = 0;
    double time, r_traj,theta_traj, time_traj, x_traj, z_traj;
    double dummy_time;
    double NL_absvol_dummmytime, meanE_absvol_dummytime;
    double NL_absvol_time,meanE_absvol_time;
    double rho_traj, P_traj, Ye_traj, v_traj, T_traj, mue_traj, mup_traj, mun_traj;
    double tmp;
    char   *tmpfname;

    FILE *fr1, *fr2, *fw1, *fw2;
    abs_dat xnu_abs_dummytime_data;
    abs_dat enu_abs_dummytime_data;
    abs_dat aenu_abs_dummytime_data;
    abs_dat xnu_abs_time_data;
    abs_dat enu_abs_time_data;
    abs_dat aenu_abs_time_data;
    bin_polated binout;
    values_at_loc xnu_trace_return;
    values_at_loc enu_trace_return;
    values_at_loc aenu_trace_return;

    fw2 = fopen(out_traj_PATH,"w");
    fr2 = fopen(in_traj_PATH, "r");
    fscanf(fr2,"%lf%lf%lf%lf%lf%lf%lf",&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp); // Reading first line of the .dat file.

    iter = 0;                   // while loop iterator.
    inrange_iter = 0;           // Count number of trajectory points within the range for a given trajdata file.
    should_open_file = false;   // File is saved only if the value of should_open_file is true. else removed.
    caught_nan = false;

    while((fscanf(fr2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &time_traj,&rho_traj,&P_traj,&Ye_traj,&v_traj,\
        &T_traj,&mue_traj,&mup_traj,&mun_traj,&r_traj,&theta_traj)) != EOF){

        iter++;
        printf("<file#: %d> %s loop iter:  %d %s\n", filecount,filler,iter,filler);
        x_traj = (r_traj*sin(theta_traj*M_PI/180.0));
        z_traj = (r_traj*cos(theta_traj*M_PI/180.0));
        time   = time_traj - trajtime_offset;

        if( ((time>=0.0)&&(time<=10.0)) && ((r_traj>=50.0)&&(r_traj<=100.0)) ){

            dummy_time = time;

            if(time<2.5){
                time = 2.5;
            }
            if(time>10.0){
                time = 10.0;
            }

            should_open_file = true;
            inrange_iter++;

            //---------------------------------------------------------------------------------------------------------

            printf("\n %s xnu %s \n\n", filler, filler);
            xnu_trace_return = Manager("xnu",time ,x_traj, z_traj);

            if(dummy_time<2.5){
                xnu_abs_dummytime_data = get_L_abs_vol(dummy_time, "xnu");
                NL_absvol_dummmytime   = xnu_abs_dummytime_data.NL_abs_vol;
                meanE_absvol_dummytime = xnu_abs_dummytime_data.meanE_abs_vol;

                xnu_abs_time_data = get_L_abs_vol(time, "xnu");
                NL_absvol_time    = xnu_abs_time_data.NL_abs_vol;
                meanE_absvol_time = xnu_abs_time_data.meanE_abs_vol;

                xnu_trace_return.n_nu    = xnu_trace_return.n_nu*(NL_absvol_dummmytime/NL_absvol_time);
                xnu_trace_return.av_e_nu = xnu_trace_return.av_e_nu*(NL_absvol_dummmytime/NL_absvol_time);
            }

            //---------------------------------------------------------------------------------------------------------

            printf("\n %s enu %s \n\n", filler, filler);
            enu_trace_return = Manager("enu",time ,x_traj, z_traj);

            if(dummy_time<2.5){
                enu_abs_dummytime_data = get_L_abs_vol(dummy_time, "enu");
                NL_absvol_dummmytime   = enu_abs_dummytime_data.NL_abs_vol;
                meanE_absvol_dummytime = enu_abs_dummytime_data.meanE_abs_vol;

                enu_abs_time_data = get_L_abs_vol(time, "enu");
                NL_absvol_time    = enu_abs_time_data.NL_abs_vol;
                meanE_absvol_time = enu_abs_time_data.meanE_abs_vol;

                enu_trace_return.n_nu    = enu_trace_return.n_nu*(NL_absvol_dummmytime/NL_absvol_time);
                enu_trace_return.av_e_nu = enu_trace_return.av_e_nu*(NL_absvol_dummmytime/NL_absvol_time);
            }

            //---------------------------------------------------------------------------------------------------------

            printf("\n %s aenu %s \n\n",  filler, filler);
            aenu_trace_return = Manager("aenu",time ,x_traj, z_traj);

            if(dummy_time<2.5){
                aenu_abs_dummytime_data = get_L_abs_vol(dummy_time, "aenu");
                NL_absvol_dummmytime    = aenu_abs_dummytime_data.NL_abs_vol;
                meanE_absvol_dummytime  = aenu_abs_dummytime_data.meanE_abs_vol;

                aenu_abs_time_data = get_L_abs_vol(time, "aenu");
                NL_absvol_time     = aenu_abs_time_data.NL_abs_vol;
                meanE_absvol_time  = aenu_abs_time_data.meanE_abs_vol;

                aenu_trace_return.n_nu    = aenu_trace_return.n_nu*(NL_absvol_dummmytime/NL_absvol_time);
                aenu_trace_return.av_e_nu = aenu_trace_return.av_e_nu*(NL_absvol_dummmytime/NL_absvol_time);
            }

            //---------------------------------------------------------------------------------------------------------

            printf("\n %s Searching in Bins %s \n\n",  filler, filler);
            binout = binterpole(dummy_time, x_traj, z_traj);
            
            //---------------------------------------------------------------------------------------------------------

            fprintf(fw2,"%le\t%le\t%le\t%le\t%le\t", dummy_time, r_traj, theta_traj, x_traj, z_traj); // col: (0-4)
            fprintf(fw2,"%le\t%le\t%le\t%le\t", rho_traj, P_traj, Ye_traj, v_traj); // col: (5-8)
            fprintf(fw2,"%le\t%le\t%le\t%le\t", T_traj, mue_traj, mup_traj, mun_traj); // col: (9-12)

            fprintf(fw2,"%le\t%le\t%le\t%le\t", binout.Nnu_rt, binout.avEnu_rt, binout.Nanu_rt, binout.avEanu_rt); // col: (13-16)

            if(enu_trace_return.n_nu == 0.0){
                enu_trace_return.n_nu = 0.0;
                enu_trace_return.av_e_nu = 0.0;
                caught_nan = true;
            }
            fprintf(fw2,"%le\t%le\t", enu_trace_return.n_nu, enu_trace_return.av_e_nu); // col: (17-18)

            if(aenu_trace_return.n_nu = 0.0){
                aenu_trace_return.n_nu = 0.0;
                aenu_trace_return.av_e_nu = 0.0;
                caught_nan = true;
            }
            fprintf(fw2,"%le\t%le\t", aenu_trace_return.n_nu, aenu_trace_return.av_e_nu); // col: (19-20)

            fprintf(fw2,"%le\t%le",   xnu_trace_return.n_nu, xnu_trace_return.av_e_nu); // col: (21-22)
            
            fprintf(fw2, "\n");

            //---------------------------------------------------------------------------------------------------------

            printf("\n<file#: %d>  %s  inrange iter: %d %s \n\n\v", filecount,filler,inrange_iter, filler);
        }
    }

    fclose(fw2);
    fclose(fr2);

    if(should_open_file == false){
        if (remove(out_traj_PATH)==0){
            printf("Removed %s\n", out_traj_PATH);
        }
    }
    
    return(caught_nan);
}
