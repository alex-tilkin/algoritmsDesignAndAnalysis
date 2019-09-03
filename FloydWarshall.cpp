//
// Created by Roie Danino on 2019-09-01.
//



#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#define INFINITY 1000000

#define ARR_SIZE 30

#define N 500

#define ROOT 0

void generateNumbers(int* arr, int size);

void printArr(const int* arr, int size);

void generateGraph(int distanceMat[N][N], int n);

void floydWarshall(int graph[N][N], int dist[N][N]);

void floydWarshallParallel(int graph[N][N], int dist[N][N], int numOfProcesses, int rank);

void printGraph(int graph[N][N]);

int main(int argc, char* argv[])
{
    int rank, numOfProcesses;
    int graph[N][N];
    int dist[N][N];
    double t0, t1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

    srand(time(NULL));


    if(rank == ROOT)
    {

        generateGraph(graph, N);
        //    printGraph(graph);
        t0 = MPI_Wtime();
        floydWarshall(graph, dist);
        t1 = MPI_Wtime();
        printf("\nTime: %lf, Shortest Paths:\n", t1 - t0);
        // printGraph(dist);
    }

    t0 =  MPI_Wtime();
    floydWarshallParallel(graph, dist, numOfProcesses, rank);
    t1 = MPI_Wtime();
    if(rank == ROOT)
    {
        printf("\nTime: %lf, Shortest Paths Parallel:\n", t1 - t0);
        //printGraph(dist);
    }


    MPI_Finalize();
    return 0;
}


void generateNumbers(int* arr, int size)
{
    int i;

    for (i = 0; i < size; i++)
    {

        arr[i] = rand() % N * N == 0 ? INFINITY : rand() % ARR_SIZE + 1;
    }
}

void generateGraph(int distanceMat[N][N], int n)
{

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        generateNumbers(distanceMat[i], n);
        distanceMat[i][i] = 0;
    }
}

void floydWarshallParallel(int graph[N][N], int dist[N][N], int numOfProcesses, int rank)
{
    int i, j, k, temp;
    int results[N][N];
    memcpy(dist, graph, sizeof(int) * N * N);

    MPI_Bcast(*dist, N * N, MPI_INT, ROOT, MPI_COMM_WORLD);


    for (k = rank; k < N; k += numOfProcesses)
    {
#pragma omp parallel for private (j)
        for (i = 0; i < N; ++i)
        {
            for (j = 0; j < N; ++j)
            {
                if(dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }


    MPI_Reduce(*dist, *results, N * N, MPI_INT, MPI_MIN, ROOT, MPI_COMM_WORLD);

    if(rank == ROOT)
        memcpy(results, dist, sizeof(int) * N * N);
}


void floydWarshall(int graph[N][N], int dist[N][N])
{
    int i, j, k;
    memcpy(dist, graph, sizeof(int) * N * N);

//#pragma omp parallel for private (i, j)
    for (k = 0; k < N; ++k)
    {
        for (i = 0; i < N; ++i)
        {
            for (j = 0; j < N; ++j)
            {
                if(dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
}

void printArr(const int* arr, int size)
{

    int i;
    for (i = 0; i < size; ++i)
    {
        if(i % ARR_SIZE == 0)
            printf("\n");
        printf("%d ", arr[i]);

    }
    printf("\n");
}

void printGraph(int graph[N][N])
{
    printf("%7s"," | ");
    for (int j = 0; j < N; ++j)
    {
        printf("%4d | ", j);
    }
    printf("\n");

    for (int j = 0; j < N - 1; ++j)
    {
        printf("- - - - - - ");
    }
    printf("\n");

    for (int i = 0; i < N; ++i)
    {
        printf("%4d | ", i);
        for (int j = 0; j < N; ++j)
        {
            if(graph[i][j] == INFINITY)
            {
                printf("%4s | ","Inf");
            }
            else
            {
                printf("%4d | ", graph[i][j]);
            }
        }
        printf("\n\n");
    }
}