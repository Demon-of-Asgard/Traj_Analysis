#include "../Libs/Libs.h"

int main(void){
    FILE *fr1, *fr2, *fw1, *fw2;

    bool caught_nan = false;

    char tmpname[300];
    char in_traj_filename[200];
    char in_traj_folder[200] = "../dd2_135_135_nu16420/";
    char in_traj_PATH[200];

    char out_traj_folder[200] = "../OutFiles/Traj_w_xnu_all_50_100/";
    char out_traj_PATH[200];

    char caughtnan_PATH[200] = "../OutFiles/caughtnan.txt";

    int filecount;

    filecount = 0;
    fw1 = fopen(caughtnan_PATH,"w");
    fr1 = fopen("../OutFiles/Traj_name_files/trajectory_names_full.txt","r");
    while(fscanf(fr1, "%s",in_traj_filename)!=EOF){

        filecount++;
        strcpy(in_traj_PATH, in_traj_folder);
        //strcat(in_traj_PATH, "trajectory.dat");
        strcat(in_traj_PATH, in_traj_filename);

        strcpy(out_traj_PATH, out_traj_folder);
        //strcat(out_traj_PATH, "trajectory.dat");
        strcat(out_traj_PATH, in_traj_filename);

        printf("in PATH: %s\n", in_traj_PATH);
        printf("out PATH: %s\n", out_traj_PATH);

        caught_nan = master(in_traj_PATH, out_traj_PATH, filecount);

        if(caught_nan == true){
            fprintf(fw1,"%s\n",in_traj_filename);
        }
    }
    fclose(fr1);
    fclose(fw1);
    return(0);
}
