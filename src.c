#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    double t1 = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    const char* filename = argv[1];

    long long PX = atoi(argv[2]), PY = atoi(argv[3]), PZ = atoi(argv[4]);
    long long NX = atoi(argv[5]), NY = atoi(argv[6]), NZ = atoi(argv[7]), NC = atoi(argv[8]);

    const char * out_filename = argv[9];

    if(argc != 10){
        perror("Argument Number is not Correct\n");
        exit(EXIT_FAILURE);
    }

    int tot_P = PX*PY*PZ;
    int tot_data = NX*NY*NZ;

    long long nx = NX/PX;
    long long ny = NY/PY;
    long long nz = NZ/PZ;

    long long leader_domain_size = NX*NY*NZ*NC;
    float * buffer = NULL;
    

    long long dom_size = nx*ny*nz*NC;
    float * domain = (float*) malloc(sizeof(float) * dom_size );


    if(rank == 0){

        MPI_File fh;

        MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);  

        MPI_Status status;

        float * leader_domain  = (float*) malloc(sizeof(float) * leader_domain_size);
        buffer =  (float*) malloc(sizeof(float) * leader_domain_size);

        MPI_Offset offset = 0;
        
        MPI_File_read_at(fh, offset, leader_domain, leader_domain_size, MPI_FLOAT, &status);

        MPI_File_close(&fh);


        long long index_of_buffer = 0;

        
        for(int p = 0; p < (PX * PY* PZ); p++ ){

            int st_x = (p%PX)*nx;
            int st_y = ((p/PX)%PY)*ny;
            int st_z = (p/(PX*PY))*nz;

            for(long long k = st_z; k<st_z + nz; k++){
                for(long long j = st_y; j < st_y + ny; j++){
                    for(long long i = st_x; i < st_x + nx; i++){

                        long long index = j * NX + i + k * NX * NY;
                        index = index*NC;

                        for(int co = 0; co < NC; co++){
                            buffer[index_of_buffer] = leader_domain[index + co];
                            index_of_buffer++;
                        }
                    }
                }
            }

        }


        if(leader_domain != NULL) free(leader_domain);


    }
    
    MPI_Scatter(buffer,dom_size,MPI_FLOAT,domain,dom_size,MPI_FLOAT,0, MPI_COMM_WORLD);

    double t2 = MPI_Wtime();
    
    if(buffer != NULL) free(buffer);

    float ****  grid = (float****) malloc(nx * sizeof(float ***));

    if(grid == NULL){
        perror("malloc failed\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0;i<nx;i++){

        grid[i] = (float***) malloc(ny * sizeof(float**));

        if(grid[i] == NULL){
            perror("malloc failed\n");
            exit(EXIT_FAILURE);
        }

        for(int j = 0; j<ny;j++){
            grid[i][j] = (float**) malloc(nz * sizeof(float*));

            if(grid[i][j] == NULL){
                perror("malloc failed\n");
                exit(EXIT_FAILURE);
            }

            for(int k = 0; k < nz; k++){
                grid[i][j][k] = (float*) malloc(NC * sizeof(float));

                if(grid[i][j][k] == NULL){
                    perror("malloc failed\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
    }


    
    long long curr = 0;
    for(int k = 0; k<nz; k++){
        for(int j = 0; j<ny; j++){
            for(int i = 0; i<nx; i++){
                for(int tc = 0; tc < NC; tc++){
                    grid[i][j][k][tc] = domain[curr];
                    curr++;
                }
            }
        }
    }


    if(domain != NULL) free(domain);

    int prcs_idx_x = rank%PX;
    int prcs_idx_y = (rank/PX) % PY;
    int prcs_idx_z = rank/(PX*PY);

    float * up_send = NULL;
    float * down_send = NULL;
    float * right_send = NULL;
    float * left_send = NULL;
    float * front_send = NULL;
    float * back_send = NULL;

    float * up_recv = NULL;
    float * down_recv = NULL;
    float * right_recv = NULL;
    float * left_recv = NULL;
    float * front_recv = NULL;
    float * back_recv = NULL;

    MPI_Request requests[12];
    int req_count = 0;

    if(prcs_idx_x > 0){
        // MPI_Irecv(right_recv);   // post Irecv from right
        // allocating in right send
        right_recv = (float *)malloc(ny*nz*NC*sizeof(float));
        MPI_Irecv(right_recv, ny * nz * NC, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        right_send = (float*) malloc(ny*nz*NC*sizeof(float));
        long long idx = 0;
        for(int k = 0; k < nz; k++){
            for(int j = 0; j< ny; j++){
                for(int tc = 0; tc < NC; tc++){
                    right_send[idx] = grid[0][j][k][tc];
                    idx++;
                }
            }
        }
    }

    if(prcs_idx_x < PX-1){
        // left recevie and left send
        // MPI_Irecv(left_recv);
        left_recv = (float *)malloc(ny * nz * NC * sizeof(float));
        MPI_Irecv(left_recv, ny * nz * NC, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        left_send = (float*) malloc(ny*nz*NC*sizeof(float));
        long long idx = 0;
        for(int k = 0; k < nz; k++){
            for(int j = 0; j< ny; j++){
                for(int tc = 0; tc < NC; tc++){
                    left_send[idx] = grid[nx-1][j][k][tc];
                    idx++;
                }
            }
        }
    }

    if(prcs_idx_y > 0){
        // MPI_Irecv(down_recv);
        down_recv = (float *)malloc(nx * nz * NC * sizeof(float));
        MPI_Irecv(down_recv, nx * nz * NC, MPI_FLOAT, rank - PX, 0, MPI_COMM_WORLD, &requests[req_count++]);
        down_send = (float *)malloc(nx * nz * NC * sizeof(float));
        long long idx = 0;
        for(int k = 0; k < nz; k++)
        {
            for(int i = 0; i < nx; i++)
            {
                for(int tc = 0; tc < NC; tc++)
                {
                    down_send[idx] = grid[i][0][k][tc];
                    idx++;
                }
            }
        }
    }

    if(prcs_idx_y < PY - 1)
    {
        // MPI_Irecv(up_recv)
        up_recv = (float *)malloc(nx * nz * NC * sizeof(float));
        MPI_Irecv(up_recv, nx * nz * NC, MPI_FLOAT, rank + PX, 0, MPI_COMM_WORLD, &requests[req_count++]);
        up_send = (float *)malloc(nx * nz * NC * sizeof(float));
        long long idx = 0;
        for(int k = 0; k < nz; k++)
        {
            for(int i = 0; i < nx; i++)
            {
                for(int tc = 0; tc < NC; tc++)
                {
                    up_send[idx] = grid[i][ny - 1][k][tc];
                    idx++;
                }
            }
        }
    }

    if(prcs_idx_z > 0)
    {
        // MPI_Irecv(back_recv)
        back_recv = (float *)malloc(nx * ny * NC * sizeof(float));
        MPI_Irecv(back_recv, nx * ny * NC, MPI_FLOAT, rank - (PX * PY), 0, MPI_COMM_WORLD, &requests[req_count++]);
        back_send = (float *)malloc(nx * ny * NC * sizeof(float));
        long long idx = 0;
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                for(int tc = 0; tc < NC; tc++)
                {
                    back_send[idx] = grid[i][j][0][tc];
                    idx++;
                }
            }
        }
    }

    if(prcs_idx_z < PZ - 1)
    {
        // MPI_Irecv(front_recv)
        front_recv = (float *)malloc(nx * ny * NC * sizeof(float));
        MPI_Irecv(front_recv, nx * ny * NC, MPI_FLOAT, rank + (PX * PY), 0, MPI_COMM_WORLD, &requests[req_count++]);

        front_send = (float *)malloc(nx * ny * NC * sizeof(float));
        long long idx = 0;
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                for(int tc = 0; tc < NC; tc++)
                {
                    front_send[idx] = grid[i][j][nz - 1][tc];
                    idx++;
                }
            }
        }
    }

    if(left_send != NULL)
    {
        MPI_Isend(left_send, ny * nz * NC, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(left_send);
    }

    if(right_send != NULL)
    {
        MPI_Isend(right_send, ny * nz * NC, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(right_send);
    }

    if(down_send != NULL)
    {
        MPI_Isend(down_send, nx * nz * NC, MPI_FLOAT, rank - PX, 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(down_send);
    }

    if(up_send != NULL)
    {
        MPI_Isend(up_send, nx * nz * NC, MPI_FLOAT, rank + PX, 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(up_send);
    }

    if(front_send != NULL)
    {
        MPI_Isend(front_send, nx * ny * NC, MPI_FLOAT, rank + (PX * PY), 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(front_send);
    }

    if(back_send != NULL)
    {
        MPI_Isend(back_send, nx * ny * NC, MPI_FLOAT, rank - (PX * PY), 0, MPI_COMM_WORLD, &requests[req_count++]);
        // free(back_send);
    }

    MPI_Status statuses[12];

    MPI_Waitall(req_count, requests, statuses);

    if(left_send == NULL) free(left_send);
    if(right_send == NULL) free(right_send);
    if(down_send == NULL) free(down_send);
    if(up_send == NULL) free(up_send);
    if(front_send == NULL) free(front_send);
    if(back_send == NULL) free(back_send);



    int * local_min_cnt = (int *)malloc(NC * sizeof(int));
    int * local_max_cnt = (int *)malloc(NC * sizeof(int));
    
    float *global_min = (float *)malloc(NC * sizeof(float));
    float *global_max = (float *)malloc(NC * sizeof(float));

    for(int i = 0; i < NC; i++)
    {
        global_max[i] = -FLT_MAX;
        global_min[i] = FLT_MAX;
        local_max_cnt[i] = 0;
        local_min_cnt[i] = 0;
    }

    for(int k = 0; k < nz; k++){
        for(int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){

                int up = j + 1;
                int down = j - 1;
                int left = i + 1;
                int right = i - 1;
                int front = k + 1;
                int back = k - 1;            

                
                for(int tc = 0; tc < NC; tc++){
                    int is_local_min = 1;
                    int is_local_max = 1;

                    if(up >= ny){
                        if(up_recv != NULL){
                            // check up
                            if(grid[i][j][k][tc] >= up_recv[(k * nx + i)*NC + tc])
                            {
                                is_local_min = 0;
                            }

                            if(grid[i][j][k][tc] <= up_recv[(k * nx + i)*NC + tc])
                            {
                                is_local_max = 0;
                            }

                        }
                    }
                    else
                    {
                        if(grid[i][j][k][tc] >= grid[i][up][k][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[i][up][k][tc])
                        {
                            is_local_max = 0;
                        }
                    }

                    if(down < 0){
                        if(down_recv != NULL){
                            // check down

                            if(grid[i][j][k][tc] >= down_recv[(k * nx + i)*NC + tc])
                            {
                                is_local_min = 0;
                            }
                            if(grid[i][j][k][tc] <= down_recv[(k * nx + i)*NC + tc])
                            {
                                is_local_max = 0;
                            }
                        }
                    }

                    else
                    {
                        if(grid[i][j][k][tc] >= grid[i][down][k][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[i][down][k][tc])
                        {
                            is_local_max = 0;
                        }
                    }

                    if(left >= nx){
                        if(left_recv != NULL){
                            // check left
                            if(grid[i][j][k][tc] >= left_recv[(k * ny + j)*NC + tc])
                            {
                                is_local_min = 0;
                            }

                            if(grid[i][j][k][tc] <= left_recv[(k * ny + j)*NC + tc])
                            {
                                is_local_max = 0;
                            }
                        }
                    }

                    else
                    {
                        if(grid[i][j][k][tc] >= grid[left][j][k][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[left][j][k][tc])
                        {
                            is_local_max = 0;
                        }
                    }


                    if(right < 0){
                        if(right_recv != NULL){
                            if(grid[i][j][k][tc] >= right_recv[(k * ny + j)*NC + tc])
                            {
                                is_local_min = 0;
                            }

                            if(grid[i][j][k][tc] <= right_recv[(k * ny + j)*NC + tc])
                            {
                                is_local_max = 0;
                            }
                        }
                    }

                    else
                    {
                        if(grid[i][j][k][tc] >= grid[right][j][k][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[right][j][k][tc])
                        {
                            is_local_max = 0;
                        }
                    }

                    if(front >= nz){
                        if(front_recv!= NULL){
                            if(grid[i][j][k][tc] >= front_recv[(j * nx + i)*NC + tc])
                            {
                                is_local_min = 0;
                            }
                            if(grid[i][j][k][tc] <= front_recv[(j * nx + i)*NC + tc])
                            {
                                is_local_max = 0;
                            }
                        }
                    }

                    else
                    {
                        if(grid[i][j][k][tc] >= grid[i][j][front][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[i][j][front][tc])
                        {
                            is_local_max = 0;
                        }
                    }


                    if(back < 0){
                        if(back_recv != NULL){
                            if(grid[i][j][k][tc] >= back_recv[(j * nx + i)*NC + tc])
                            {
                                is_local_min = 0;
                            }
                            if(grid[i][j][k][tc] <= back_recv[(j * nx + i)*NC + tc])
                            {
                                is_local_max = 0;
                            }
                        }
                    }

                    else
                    {
                        if(grid[i][j][k][tc] >= grid[i][j][back][tc])
                        {
                            is_local_min = 0;
                        }

                        if(grid[i][j][k][tc] <= grid[i][j][back][tc])
                        {
                            is_local_max = 0;
                        }
                    }


                    if(is_local_max){
                        local_max_cnt[tc]++;
                        if(grid[i][j][k][tc] > global_max[tc]){
                            global_max[tc] = grid[i][j][k][tc];
                        }
                    }

                    if(is_local_min){
                        local_min_cnt[tc]++;
                        if(grid[i][j][k][tc] < global_min[tc]){
                            global_min[tc] = grid[i][j][k][tc];
                        }
                    }

                }


            }
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                if(grid[i][j][k] != NULL) free(grid[i][j][k]);  
            }
            if(grid[i][j] != NULL) free(grid[i][j]);        
        }
        if(grid[i] != NULL) free(grid[i]);                
    }

    if(grid != NULL) free(grid);  
    
    
    free(right_recv);
    free(left_recv);

    free(up_recv);
    free(down_recv);

    free(front_recv);
    free(back_recv);



    int * total_min_cnt = NULL;
    int * total_max_cnt = NULL;
    float *final_global_min = NULL;
    float *final_global_max = NULL;

    if(rank == 0){
        total_min_cnt = (int *)malloc(NC * sizeof(int));
        total_max_cnt = (int *)malloc(NC * sizeof(int));
        final_global_min = (float *)malloc(NC * sizeof(float));
        final_global_max = (float *)malloc(NC * sizeof(float));
    }

    
    MPI_Reduce(local_min_cnt,total_min_cnt,NC,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(local_max_cnt,total_max_cnt,NC,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    

    MPI_Reduce(global_min,final_global_min,NC,MPI_FLOAT,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(global_max,final_global_max,NC,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);


   
    double t3 = MPI_Wtime();

    double computation_time = t3 - t2;
    double communication_time = t2 - t1;
    double total_time = t3 - t1;

    double final_computation_time;
    double final_communication_time;
    double final_total_time;

    MPI_Reduce(&computation_time,&final_computation_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&communication_time,&final_communication_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&total_time,&final_total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    

    if(rank == 0){


        FILE *fp = fopen(out_filename, "w");
        if (fp == NULL) {
            perror("Error opening file");
            return EXIT_FAILURE;
        }


        for(int i = 0; i<NC; i++){
            fprintf(fp,"(%d, %d)",total_min_cnt[i],total_max_cnt[i]);
            if(i != NC-1) fprintf(fp,",");
        }

        fprintf(fp,"\n");

        for(int i = 0; i<NC; i++){
            fprintf(fp,"(%.4f, %.4f)",final_global_min[i],final_global_max[i]);
            if(i != NC-1) fprintf(fp,",");
        }

        fprintf(fp,"\n%.4f %.4f %.4f",final_communication_time,final_computation_time,final_total_time);

    }
    
    
    if(local_max_cnt != NULL) free(local_max_cnt);
    if(local_min_cnt != NULL) free(local_min_cnt);

    if(global_max != NULL) free(global_max);
    if(global_min != NULL) free(global_min);
    
    if(final_global_max != NULL) free(final_global_max);
    if(final_global_min != NULL) free(final_global_min);

    if(total_max_cnt != NULL) free(total_max_cnt);
    if(total_min_cnt != NULL) free(total_min_cnt);


    MPI_Finalize();
    return 0;
}
