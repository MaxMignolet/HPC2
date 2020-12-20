#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MAX_SIZE 100
#define MAX_ITERATION 40 // for the conjugate gradient

void saveToFile(double** u, double** v, int rank, int world_size, int N,
                int n_y, int timeStep);

void getBorders(double** x, double* top_border, double* bottom_border, int rank,
                int rank_below, int rank_above, int n_y, int N);

void implicit_multiByA(double** x, double** y, double C, int rank,
                       int rank_below, int rank_above, int n_y, int N);

void implicit_matrix_add(double** a, double** b, double** c, int n_y, int N);

void implicit_matrix_sub(double** a, double** b, double** c, int n_y, int N);

double implicit_martix_innerProd(double** a, double** b, int n_y, int N);

// void free3Dtab(double ***tab, int x, int y);

void implicit_conjugate_gradient(double** x, double** b, double C, int n_y,
                                 int N, double r_treshold, double** r_i,
                                 double** p_i, double** a_p_i, double** a_x0,
                                 int rank, int rank_below, int rank_above);

int main(int argc, char** argv){
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0){fprintf(stdout, "Number of processes: %d\n", world_size);}
    if(rank == 0){fprintf(stdout, "Number of threads/process %d\n", omp_get_max_threads());}

    // starts timer
    double start_time, end_time, start_time_comp, end_time_comp;
    if(rank == 0){start_time = MPI_Wtime();}

    //get parameters from file param.txt
    FILE* fp;
    fp = fopen(argv[1], "r");
    int resolution_scheme = atoi(argv[2]);

    char tmp[MAX_SIZE] = "";
    int nb_param = 9;
    double param[nb_param-1];
    int i, j, l, t, eta, xi; // variables for for loops
    // general rule: i, j -> global coords in the main grid (x, y dir resp.)
    //               eta, xi -> local coords in the local grid (x, y dir resp.)
    //               t -> time step
    //               s -> frame number
    //               l, m -> for the other loops
    //               lambda -> chunk index along y
    // this general rule may not hold in subroutines

    if (fp == NULL){
        fprintf(stderr, "Can't open the file\n");
        return 1;
    }
    else{
        for(l = 0; l < nb_param-1; l++){
            fgets(tmp, MAX_SIZE, fp);
            param[l] = atof(tmp);
        }
    }

    double dt = param[0];
    double dx = param[1];
    double T_max = param[2];
    double D_u = param[3];
    double D_v = param[4];
    double f = param[5];
    double k = param[6];
    double r_treshold = param[7];

    fgets(tmp, MAX_SIZE, fp);
    fclose(fp);
    unsigned int S = atoi(tmp);

    // ranks of chunks above and below
    int rank_above, rank_below;
    if(rank == 0){
        rank_below = world_size - 1;
    }else{
        rank_below = rank - 1;
    }
    if(rank == world_size - 1){
        rank_above = 0;
    }else{
        rank_above = rank + 1;
    }

    // total number of elements along a direction (total grid)
    unsigned int N = (unsigned int)floor(1/dx)+1;
    if(rank == 0){fprintf(stdout, "N = %d\n", N);}
    // total number of time steps
    unsigned int T = (unsigned int)floor(T_max/dt)+1;
    if(rank == 0){fprintf(stdout, "T = %u\n", T);}

    // number of elements along the y direction IN a chunk
    unsigned int n_y = N / world_size;
    if(rank == world_size - 1){
        n_y = N - n_y * (world_size - 1);
    }

    // allocate sub-grids
    // first index = local coordinate along y
    // second index = local coordinate along x
    double **u = NULL;
    double **v = NULL;

    u = malloc(sizeof(double*) * n_y);
    v = malloc(sizeof(double*) * n_y);
    if(u == NULL || v == NULL){
        fprintf(stderr, "Allocation error, type 1\n");
        return 1;
    }

    for(xi = 0; xi < n_y; xi++){
        u[xi] = malloc(sizeof(double) * N);
        v[xi] = malloc(sizeof(double) * N);
        if(u[xi] == NULL || u[xi] == NULL){
            fprintf(stderr, "Allocation error, type 2\n");
            return 1;
        }
    }

    //initial conditions
    for(xi = 0; xi < n_y; xi++){
        // conversion from local coords to global
        j = xi + rank * (N / world_size);
        for(i = 0; i < N; i++){
            u[xi][i] = 1;
            // circle of radius 0.05 at the centre of the square domain
            if(sqrt(pow((double)i/(double)N-0.5, 2) 
                    + pow((double)j/(double)N-0.5, 2)) <= 0.05){
                v[xi][i] = 1;
            }
            else{
                v[xi][i] = 0;
            }
        }
    }

    if(S != 0){
        t = 0;
        saveToFile(u, v, rank, world_size, N, n_y, t);
    }

    if(rank == 0){start_time_comp = MPI_Wtime();}

    if(resolution_scheme == 0){ // explicit euler scheme
        // temporary variables for solving the problem
        // u_prev -> solution for u at time t-1
        // u_new  -> solution for u at time t
        double** u_prev = NULL;
        double** u_new = NULL;
        double** v_prev = NULL;
        double** v_new = NULL;

        u_prev = u;
        u_new = malloc(sizeof(double*) * n_y);
        v_prev = v;
        v_new = malloc(sizeof(double*) * n_y);

        for(xi = 0; xi < n_y; xi++){
            u_new[xi] = malloc(sizeof(double) * N);
            v_new[xi] = malloc(sizeof(double) * N);
        }

        // just to keep a reference to the allocated memory so we can free it
        // later, see lines 426 - 433 (end of explicit scheme when freeing memory)
        double** ref_2_u_2_free = u_new;
        double** ref_2_v_2_free = v_new;

        double* top_border_u = NULL;
        double* bottom_border_u = NULL;
        double* top_border_v = NULL;
        double* bottom_border_v = NULL;
        top_border_u = malloc(sizeof(double) * N);
        bottom_border_u = malloc(sizeof(double) * N);
        top_border_v = malloc(sizeof(double) * N);
        bottom_border_v = malloc(sizeof(double) * N);

        // some temporary variables
        double u_middle, u_above, u_below, u_left, u_right;
        double v_middle, v_above, v_below, v_left, v_right;
        double C_u = D_u / pow(dx, 2);
        double C_v = D_v / pow(dx, 2);
        double** temp;

        for(t = 1; t <= T; t++){
            // computing inner of chunk
            #pragma omp parallel for collapse(2)\
                    private(u_middle, u_below, u_above, u_left, u_right,\
                            v_middle, v_below, v_above, v_left, v_right)
            for(xi = 1; xi < n_y - 1; xi++){
                for(eta = 1; eta < N - 1; eta++){
                    u_middle = u_prev[xi][eta];
                    u_below = u_prev[xi-1][eta];
                    u_above = u_prev[xi+1][eta];
                    u_left = u_prev[xi][eta-1];
                    u_right = u_prev[xi][eta+1];

                    v_middle = v_prev[xi][eta];
                    v_below = v_prev[xi-1][eta];
                    v_above = v_prev[xi+1][eta];
                    v_left = v_prev[xi][eta-1];
                    v_right = v_prev[xi][eta+1];

                    u_new[xi][eta] = u_middle + dt
                        * (C_u * (u_above + u_right - 4 * u_middle + u_below
                                  + u_left)
                            - u_middle * pow(v_middle, 2)
                            + f * (1 - u_middle));
                    v_new[xi][eta] = v_middle + dt 
                        * (C_v * (v_above + v_right - 4 * v_middle + v_below
                                  + v_left)
                            + u_middle * pow(v_middle, 2)
                            - (f + k) * v_middle);
                }
            }

            // left and right borders
            #pragma omp parallel for\
                    private(u_middle, u_below, u_above, u_left, u_right,\
                            v_middle, v_below, v_above, v_left, v_right,\
                            eta)
            for(xi = 1; xi < n_y - 1; xi++){
                // left border
                eta = 0;

                u_middle = u_prev[xi][eta];
                u_below = u_prev[xi-1][eta];
                u_above = u_prev[xi+1][eta];
                u_left = u_prev[xi][n_y-1]; //
                u_right = u_prev[xi][eta+1];

                v_middle = v_prev[xi][eta];
                v_below = v_prev[xi-1][eta];
                v_above = v_prev[xi+1][eta];
                v_left = v_prev[xi][n_y-1]; //
                v_right = v_prev[xi][eta+1];

                u_new[xi][eta] = u_middle + dt
                    * (C_u * (u_above + u_right - 4 * u_middle + u_below
                              + u_left)
                        - u_middle * pow(v_middle, 2)
                        + f * (1 - u_middle));
                v_new[xi][eta] = v_middle + dt 
                    * (C_v * (v_above + v_right - 4 * v_middle + v_below
                              + v_left)
                        + u_middle * pow(v_middle, 2)
                        - (f + k) * v_middle);

                // right border
                eta = N - 1;
                
                u_middle = u_prev[xi][eta];
                u_below = u_prev[xi-1][eta];
                u_above = u_prev[xi+1][eta];
                u_left = u_prev[xi][eta-1];
                u_right = u_prev[xi][0]; //

                v_middle = v_prev[xi][eta];
                v_below = v_prev[xi-1][eta];
                v_above = v_prev[xi+1][eta];
                v_left = v_prev[xi][eta-1];
                v_right = v_prev[xi][0]; //

                u_new[xi][eta] = u_middle + dt
                    * (C_u * (u_above + u_right - 4 * u_middle + u_below
                              + u_left)
                        - u_middle * pow(v_middle, 2)
                        + f * (1 - u_middle));
                v_new[xi][eta] = v_middle + dt 
                    * (C_v * (v_above + v_right - 4 * v_middle + v_below
                              + v_left)
                        + u_middle * pow(v_middle, 2)
                        - (f + k) * v_middle);
            }

            // we have to retrieve the data at the top/bottom borders
            getBorders(u_prev, top_border_u, bottom_border_u, rank, rank_below, 
                       rank_above, n_y, N);
            getBorders(v_prev, top_border_v, bottom_border_v, rank, rank_below,
                       rank_above, n_y, N);

            // bottom and top borders
            #pragma omp parallel for\
                    private(u_middle, u_below, u_above, u_left, u_right,\
                            v_middle, v_below, v_above, v_left, v_right,\
                            xi)
            for(eta = 1; eta < N - 1; eta++){
                // bottom border
                xi = 0;

                u_middle = u_prev[xi][eta];
                u_below = bottom_border_u[eta]; //
                u_above = u_prev[xi+1][eta];
                u_left = u_prev[xi][eta-1];
                u_right = u_prev[xi][eta+1];

                v_middle = v_prev[xi][eta];
                v_below = bottom_border_v[eta]; //
                v_above = v_prev[xi+1][eta];
                v_left = v_prev[xi][eta-1];
                v_right = v_prev[xi][eta+1];

                u_new[xi][eta] = u_middle + dt
                    * (C_u * (u_above + u_right - 4 * u_middle + u_below
                              + u_left)
                        - u_middle * pow(v_middle, 2)
                        + f * (1 - u_middle));
                v_new[xi][eta] = v_middle + dt 
                    * (C_v * (v_above + v_right - 4 * v_middle + v_below
                              + v_left)
                        + u_middle * pow(v_middle, 2)
                        - (f + k) * v_middle);

                // top border
                xi = n_y - 1;
                
                u_middle = u_prev[xi][eta];
                u_below = u_prev[xi-1][eta];
                u_above = top_border_u[eta]; //
                u_left = u_prev[xi][eta-1];
                u_right = u_prev[xi][eta+1];

                v_middle = v_prev[xi][eta];
                v_below = v_prev[xi-1][eta];
                v_above = top_border_v[eta]; //
                v_left = v_prev[xi][eta-1];
                v_right = v_prev[xi][eta+1];

                u_new[xi][eta] = u_middle + dt
                    * (C_u * (u_above + u_right - 4 * u_middle + u_below
                              + u_left)
                        - u_middle * pow(v_middle, 2)
                        + f * (1 - u_middle));
                v_new[xi][eta] = v_middle + dt 
                    * (C_v * (v_above + v_right - 4 * v_middle + v_below
                              + v_left)
                        + u_middle * pow(v_middle, 2)
                        - (f + k) * v_middle);
            }

            // corners now
            // bottom left
            xi = 0;
            eta = 0;

            u_middle = u_prev[xi][eta];
            u_below = bottom_border_u[eta]; //
            u_above = u_prev[xi+1][eta];
            u_left = u_prev[xi][N-1]; //
            u_right = u_prev[xi][eta+1];

            v_middle = v_prev[xi][eta];
            v_below = bottom_border_v[eta]; //
            v_above = v_prev[xi+1][eta];
            v_left = v_prev[xi][N-1]; //
            v_right = v_prev[xi][eta+1];

            u_new[xi][eta] = u_middle + dt
                * (C_u * (u_above + u_right - 4 * u_middle + u_below
                          + u_left)
                    - u_middle * pow(v_middle, 2)
                    + f * (1 - u_middle));
            v_new[xi][eta] = v_middle + dt 
                * (C_v * (v_above + v_right - 4 * v_middle + v_below
                          + v_left)
                    + u_middle * pow(v_middle, 2)
                    - (f + k) * v_middle);

            // bottom right
            xi = 0;
            eta = N - 1;
            u_middle = u_prev[xi][eta];
            u_below = bottom_border_u[eta]; //
            u_above = u_prev[xi+1][eta];
            u_left = u_prev[xi][eta-1];
            u_right = u_prev[xi][0]; //

            v_middle = v_prev[xi][eta];
            v_below = bottom_border_v[eta]; //
            v_above = v_prev[xi+1][eta];
            v_left = v_prev[xi][eta-1];
            v_right = v_prev[xi][0]; //

            u_new[xi][eta] = u_middle + dt
                * (C_u * (u_above + u_right - 4 * u_middle + u_below
                          + u_left)
                    - u_middle * pow(v_middle, 2)
                    + f * (1 - u_middle));
            v_new[xi][eta] = v_middle + dt 
                * (C_v * (v_above + v_right - 4 * v_middle + v_below
                          + v_left)
                    + u_middle * pow(v_middle, 2)
                    - (f + k) * v_middle);

            // top left
            xi = n_y - 1;
            eta = 0;
            u_middle = u_prev[xi][eta];
            u_below = u_prev[xi-1][eta];
            u_above = top_border_u[eta]; //
            u_left = u_prev[xi][N-1]; //
            u_right = u_prev[xi][eta+1];

            v_middle = v_prev[xi][eta];
            v_below = v_prev[xi-1][eta];
            v_above = top_border_v[eta]; //
            v_left = v_prev[xi][N-1]; //
            v_right = v_prev[xi][eta+1];

            u_new[xi][eta] = u_middle + dt
                * (C_u * (u_above + u_right - 4 * u_middle + u_below
                          + u_left)
                    - u_middle * pow(v_middle, 2)
                    + f * (1 - u_middle));
            v_new[xi][eta] = v_middle + dt 
                * (C_v * (v_above + v_right - 4 * v_middle + v_below
                          + v_left)
                    + u_middle * pow(v_middle, 2)
                    - (f + k) * v_middle);

            // top right
            xi = n_y - 1;
            eta = N - 1;
            u_middle = u_prev[xi][eta];
            u_below = u_prev[xi-1][eta];
            u_above = top_border_u[eta]; //
            u_left = u_prev[xi][eta-1];
            u_right = u_prev[xi][0]; //

            v_middle = v_prev[xi][eta];
            v_below = v_prev[xi-1][eta];
            v_above = top_border_v[eta]; //
            v_left = v_prev[xi][eta-1];
            v_right = v_prev[xi][0]; //

            u_new[xi][eta] = u_middle + dt
                * (C_u * (u_above + u_right - 4 * u_middle + u_below
                          + u_left)
                    - u_middle * pow(v_middle, 2)
                    + f * (1 - u_middle));
            v_new[xi][eta] = v_middle + dt 
                * (C_v * (v_above + v_right - 4 * v_middle + v_below
                          + v_left)
                    + u_middle * pow(v_middle, 2)
                    - (f + k) * v_middle);

            // save u and v if needed
            if(S != 0 && t % S == 0){
                // if(rank == 0){fprintf(stderr, "Saving to file, timeStep = %d,...\n", t);}
                saveToFile(u_new, v_new, rank, world_size, N, n_y, t);
                // if(rank == 0){fprintf(stderr, "Saved\n");}
            }
            // the new solution becomes the old one
            // memory allocated to old solution gets recycled
            temp = u_prev;
            u_prev = u_new;
            u_new = temp;

            temp = v_prev;
            v_prev = v_new;
            v_new = temp;
        }
        // let us free what we don't need anymore:
        for(xi = 0; xi < n_y; xi++){
            free(ref_2_u_2_free[xi]);
            free(ref_2_v_2_free[xi]);
        }
        free(ref_2_u_2_free);
        free(ref_2_v_2_free);
    }
    else if(resolution_scheme == 1){ // semi-implicit scheme
        // TODO: add a description for each variable, and explain what's going on here
        // some temporary variables
        double C_u = dt * D_u / pow(dx, 2);
        double C_v = dt * D_v / pow(dx, 2);

        double** u_i = u; // it's just a change of name to agree with the
        double** v_i = v; // notations in the pseudo code of the statement

        // temporary variables for the conjugate gradient / matrices
        double** b_u = NULL;
        double** b_v = NULL;
        double** r_i = NULL;
        double** p_i = NULL;
        double** a_x0 = NULL;
        double** a_p_i = NULL;

        b_u = malloc(sizeof(double*) * n_y);
        b_v = malloc(sizeof(double*) * n_y);
        r_i = malloc(sizeof(double*) * n_y);
        p_i = malloc(sizeof(double*) * n_y);
        a_x0 = malloc(sizeof(double*) * n_y);
        a_p_i = malloc(sizeof(double*) * n_y);

        for(xi = 0; xi < n_y; xi++){
            b_u[xi] = malloc(sizeof(double) * N);
            b_v[xi] = malloc(sizeof(double) * N);
            r_i[xi] = malloc(sizeof(double) * N);
            p_i[xi] = malloc(sizeof(double) * N);
            a_x0[xi] = malloc(sizeof(double) * N);
            a_p_i[xi] = malloc(sizeof(double) * N);
        }

        for(t = 1; t <= T; t++){
            // first let us compute the right-hand side b
            #pragma omp parallel for collapse(2)
            for(xi = 0; xi < n_y; xi++){
                for(eta = 0; eta < N; eta++){
                    b_u[xi][eta] = u_i[xi][eta]
                        + dt * (-u_i[xi][eta] * pow(v_i[xi][eta], 2) 
                                + f * (1 - u_i[xi][eta]));
                    b_v[xi][eta] = v_i[xi][eta]
                        + dt * (u_i[xi][eta] * pow(v_i[xi][eta], 2) 
                                - (f + k) * v_i[xi][eta]);
                }
            }

            // gradient conjugate: variables:
                // u_i, b_u, C_u, a_x0, r_i, p_i, a_p_i, rank, rank_below,
                // rank_above, n_y, N
            implicit_conjugate_gradient(u_i, b_u, C_u, n_y, N, r_treshold, r_i,
                                        p_i, a_p_i, a_x0, rank, rank_below,
                                        rank_above);
            implicit_conjugate_gradient(v_i, b_v, C_v, n_y, N, r_treshold, r_i,
                                        p_i, a_p_i, a_x0, rank, rank_below,
                                        rank_above);
            
            // save u and v if needed
            if(S != 0 && t % S == 0){
                saveToFile(u_i, v_i, rank, world_size, N, n_y, t);
            }
        }

        // let us free what we don't need anymore:
        // to be completed...
        for(xi = 0; xi < n_y; xi++){
            free(b_u[xi]);
            free(b_v[xi]);
            free(r_i[xi]);
            free(p_i[xi]);
            free(a_x0[xi]);
            free(a_p_i[xi]);
        }
        free(b_u);
        free(b_v);
        free(r_i);
        free(p_i);
        free(a_x0);
        free(a_p_i);
    }

    if(rank == 0){
        end_time_comp = MPI_Wtime();
        fprintf(stdout, "Elapsed wall-time (computational part): %.4lf seconds\n",
                end_time_comp - start_time_comp);
    }

    //FREE, TODO
    for(xi = 0; xi < n_y; xi++){
        free(u[xi]);
        free(v[xi]);
    }
    free(u);
    free(v);

    if(rank == 0){
        end_time = MPI_Wtime();
        fprintf(stdout, "Elapsed wall-time (total): %.4lf seconds\n",
                end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}

void saveToFile(double** u, double** v, int rank, int world_size, int N,
                int n_y, int timeStep){
    int i, j, xi, lambda;
    MPI_Status status;
    // gather all the solutions towards node with rank 0
    // gather only one time step every S time step, there's no need to gather 
    // the data that will not be saved
    double **u_to_be_saved = NULL;
    double **v_to_be_saved = NULL;

    // buffers (used only by the root)
    double* recvbuf_u = NULL;
    double* sendbuf_u = NULL;
    double* recvbuf_v = NULL;
    double* sendbuf_v = NULL;
    int root = 0; // let us avoid magic numbers

    if(rank == 0){ // some init before gathering the data
        // allocate grids for the whole domain
        u_to_be_saved = malloc(sizeof(double**) * N);
        v_to_be_saved = malloc(sizeof(double**) * N);
        if(u_to_be_saved == NULL || v_to_be_saved == NULL){
            fprintf(stderr, "Allocation error, type 1 - saveToFile\n");
            return;
        }

        for(i = 0; i < N; i++){
            u_to_be_saved[i] = malloc(sizeof(double) * N);
            v_to_be_saved[i] = malloc(sizeof(double) * N);
            if(u_to_be_saved[i] == NULL || v_to_be_saved[i] == NULL){
                fprintf(stderr, "Allocation error, type 2 - saveToFile\n");
                return;
            }
        }
    }

    for(j = 0; j < N; j++){
        if(rank != world_size - 1){
            lambda = j / n_y; // index of the chunk we retrieve data from
            if(lambda == world_size){
                lambda = world_size - 1;
            }
            xi = j - lambda * n_y; // local coord along y
        }else{
            lambda = j / (N / world_size);
            if(lambda == world_size){
                lambda = world_size - 1;
            }
            xi = j - lambda * (N / world_size);
        }

        if(rank == root){
            recvbuf_u = u_to_be_saved[j];
            recvbuf_v = v_to_be_saved[j];
            if(lambda != root){
                MPI_Recv(recvbuf_u, N, MPI_DOUBLE, lambda, 1,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(recvbuf_v, N, MPI_DOUBLE, lambda, 1,
                         MPI_COMM_WORLD, &status);
            }else{
                for(i = 0; i < N; i++){
                    recvbuf_u[i] = u[xi][i];
                    recvbuf_v[i] = v[xi][i];
                }
            }
        }else{
            if(rank == lambda){
                sendbuf_u = u[xi];
                sendbuf_v = v[xi];
                MPI_Send(sendbuf_u, N, MPI_DOUBLE, root, 1,
                         MPI_COMM_WORLD);
                MPI_Send(sendbuf_v, N, MPI_DOUBLE, root, 1,
                         MPI_COMM_WORLD);
            }
        }
    }

    if(rank == 0){ // save the solutions in binary files
        // a subdirectory named "results" must exist beforehand
        FILE* fp_u;
        FILE* fp_v;
        char t_char[10];
        sprintf(t_char, "%u", timeStep);

        char filename_u[25] = "results/u_";
        strcat(filename_u, t_char);
        strcat(filename_u, ".dat");
        char filename_v[25] = "results/v_";
        strcat(filename_v, t_char);
        strcat(filename_v, ".dat");

        fp_u = fopen(filename_u, "wb");
        if(fp_u == NULL){
            fprintf(stderr, "Error while writing data - saveToFile\n");
            return;
        }
        fp_v = fopen(filename_v, "wb");
        if(fp_u == NULL){
            fprintf(stderr, "Error while writing data - saveToFile\n");
            return;
        }
        fwrite(&N, sizeof(N), 1, fp_u);
        fwrite(&N, sizeof(N), 1, fp_v);
        for(i = 0; i < N; i++){
            fwrite(u_to_be_saved[i], sizeof(double), N, fp_u);
            fwrite(v_to_be_saved[i], sizeof(double), N, fp_v);
        }
        fclose(fp_u);
        fclose(fp_v);
    }

    return;
}

void getBorders(double** x, double* top_border, double* bottom_border, int rank,
                int rank_below, int rank_above, int n_y, int N){
    /* retrieve top and bottom borders from neighbours
     * u and v: (n_y) * (N) martices
     */
    MPI_Status status;
    
    // processes with an odd rank, send their top border, then receive
    // processes with an even rank, receive their bottom border, then send
    if(rank != rank_below){ // if(world_size != 1)
        if(rank % 2){ // odd rank
            MPI_Send(x[n_y-1], N, MPI_DOUBLE, rank_above, 1, MPI_COMM_WORLD);
            MPI_Recv(top_border, N, MPI_DOUBLE, rank_above, 1, MPI_COMM_WORLD,
                     &status);
        }else{ // even rank
            MPI_Recv(bottom_border, N, MPI_DOUBLE, rank_below, 1, MPI_COMM_WORLD,
                     &status);
            MPI_Send(x[0], N, MPI_DOUBLE, rank_below, 1, MPI_COMM_WORLD);
        }
        // afterwards same in the other way around
        if(rank % 2){ // odd rank
            MPI_Recv(bottom_border, N, MPI_DOUBLE, rank_below, 1, MPI_COMM_WORLD,
                     &status);
            MPI_Send(x[0], N, MPI_DOUBLE, rank_below, 1, MPI_COMM_WORLD);
        }
        else{ // even rank
            MPI_Send(x[n_y-1], N, MPI_DOUBLE, rank_above, 1, MPI_COMM_WORLD);
            MPI_Recv(top_border, N, MPI_DOUBLE, rank_above, 1, MPI_COMM_WORLD,
                     &status);
        }
    }else{ // if there is only one process
        #pragma omp parallel for
        for(int eta = 0; eta < N; eta++){
            top_border[eta] = x[0][eta];
            bottom_border[eta] = x[n_y-1][eta];
        }
    }
    return;
}

void implicit_multiByA(double** x, double** y, double C, int rank,
                       int rank_below, int rank_above, int n_y, int N){
    /* computes y = A x
     * y: n_y * N matrix
     * x: n_y * N matrix
     * A: 4-th order tensor representing the linear system to be solved for
     *     the semi-implicit scheme
     */
    int i, j;
    double x_middle, x_above, x_below, x_left, x_right;
    double* bottom_border = NULL;
    double* top_border = NULL;

    #pragma omp parallel for collapse(2)\
            private(x_middle, x_below, x_above, x_left, x_right)
    for(i = 1; i < n_y-1; i++){ // inner part of x and y (borders are excluded)
        for(j = 1; j < N-1; j++){
            x_middle = x[i][j];
            x_above = x[i+1][j];
            x_below = x[i-1][j];
            x_right = x[i][j+1];
            x_left = x[i][j-1];
            y[i][j] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
        }
    }
    // right and left borders
    #pragma omp parallel for\
            private(x_middle, x_below, x_above, x_left, x_right, j)
    for(i = 1; i < n_y-1; i++){ // use openMP here ?
        // left border
        j = 0;
        x_middle = x[i][j];
        x_above = x[i+1][j];
        x_below = x[i-1][j];
        x_right = x[i][j+1];
        x_left = x[i][n_y-1];
        y[i][j] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
        // right border
        j = N-1;
        x_middle = x[i][j];
        x_above = x[i+1][j];
        x_below = x[i-1][j];
        x_right = x[i][0];
        x_left = x[i][j-1];
        y[i][j] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
    }

    // top and bottom parts
    // retrieve borders from top and bottom neighbours
    top_border = malloc(sizeof(double) * N);
    bottom_border = malloc(sizeof(double) * N);
    getBorders(x, top_border, bottom_border, rank, rank_below, rank_above, n_y, N);

    // we can now compute the top and bottom parts of x
    #pragma omp parallel for\
            private(x_middle, x_below, x_above, x_left, x_right, i)
    for(j = 1; j < N-1; j++){ // use openMP here
        // top part
        i = n_y-1;
        x_middle = x[i][j];
        x_above = top_border[j];
        x_below = x[i-1][j];
        x_right = x[i][j+1];
        x_left = x[i][j-1];
        y[i][j] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
        // bottom part
        i = 0;
        x_middle = x[i][j];
        x_above = x[i+1][j];
        x_below = bottom_border[j];
        x_right = x[i][j+1];
        x_left = x[i][j-1];
        y[i][j] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
    }

    // let us compute the corners while we are at it
    // bottom left corner (0, 0)
    x_middle = x[0][0];
    x_above = x[1][0];
    x_below = bottom_border[0];
    x_right = x[0][1];
    x_left = x[0][N-1];
    y[0][0] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
    // bottom right corner (0, N-1)
    x_middle = x[0][N-1];
    x_above = x[1][N-1];
    x_below = bottom_border[N-1];
    x_right = x[0][0];
    x_left = x[0][N-2];
    y[0][N-1] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
    // top left corner (n_y-1, 0)
    x_middle = x[n_y-1][0];
    x_above = top_border[0];
    x_below = x[n_y-2][0];
    x_right = x[n_y-1][1];
    x_left = x[n_y-1][N-1];
    y[n_y-1][0] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);
    // top right corner (n_y-1, N-1)
    x_middle = x[n_y-1][N-1];
    x_above = top_border[N-1];
    x_below = x[n_y-2][N-1];
    x_right = x[n_y-1][0];
    x_left = x[n_y-1][N-2];
    y[n_y-1][N-1] = x_middle - C * (x_left + x_right - 4*x_middle + x_below + x_above);

    if(rank != rank_below){ // clear memory if it was allocated
        free(top_border);
        free(bottom_border);
    }

    return;
}

void implicit_matrix_add(double** a, double** b, double** c, int n_y, int N){
    // c = a + b
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n_y; i++){ // use openMP here
        for(int j = 0; j < N; j++){
            c[i][j] = a[i][j] + b[i][j];
        }
    }
}

void implicit_matrix_sub(double** a, double** b, double** c, int n_y, int N){
    // c = a - b
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n_y; i++){ // use openMP here
        for(int j = 0; j < N; j++){
            c[i][j] = a[i][j] - b[i][j];
        }
    }
}

double implicit_martix_innerProd(double** a, double** b, int n_y, int N){
    // ret = sum_ij a_ij * b_ij
    double local_innerProd = 0;
    double global_innerProd = 0;
    #pragma omp parallel for collapse(2) reduction(+:local_innerProd)
    for(int i = 0; i < n_y; i++){ // use openMP here
        for(int j = 0; j < N; j++){
            local_innerProd += a[i][j] * b[i][j];
        }
    }
    MPI_Allreduce(&local_innerProd, &global_innerProd, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    return global_innerProd;
}

void implicit_conjugate_gradient(double** x, double** b, double C, int n_y,
                                 int N, double r_treshold, double** r_i,
                                 double** p_i, double** a_p_i, double** a_x0,
                                 int rank, int rank_below, int rank_above){
    int eta, xi;
    double alpha, beta;
    double r_norm_0, r_norm_i, r_norm_i_1, pAp_i;
    int iteration_Nb;
    implicit_multiByA(x, a_x0, C, rank, rank_below, rank_above,
                      n_y, N);
    implicit_matrix_sub(b, a_x0, r_i, n_y, N);
    #pragma omp parallel for collapse(2)
    for(xi = 0; xi < n_y; xi++){
        for(eta = 0; eta < N; eta++){
            p_i[xi][eta] = r_i[xi][eta];
        }
    }
    r_norm_0 = implicit_martix_innerProd(r_i, r_i, n_y, N);
    r_norm_i = r_norm_0;
    iteration_Nb = 0;
    while(iteration_Nb < MAX_ITERATION){
        implicit_multiByA(p_i, a_p_i, C, rank, rank_below, rank_above, n_y, N);
        pAp_i = implicit_martix_innerProd(p_i, a_p_i, n_y, N);
        alpha = r_norm_i / pAp_i;
        
        #pragma omp parallel for collapse(2)
        for(xi = 0; xi < n_y; xi++){
            for(eta = 0; eta < N; eta++){
                x[xi][eta] += alpha * p_i[xi][eta];
                r_i[xi][eta] -= alpha * a_p_i[xi][eta];
            }
        }
        
        r_norm_i_1 = implicit_martix_innerProd(r_i, r_i, n_y, N);

        if(r_norm_i_1 < pow(r_norm_0, 2) * r_treshold){
            break;
        }

        beta = r_norm_i_1 / r_norm_i;
        
        #pragma omp parallel for collapse(2)
        for(xi = 0; xi < n_y; xi++){
            for(eta = 0; eta < N; eta++){
                p_i[xi][eta] = r_i[xi][eta] + beta * p_i[xi][eta];
            }
        }
        r_norm_i = r_norm_i_1;
        iteration_Nb++;
    }
    return;
}
