#include<mpi.h>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include <algorithm>

long int startPos;
long int endPos;
long int partSize;
long int rightPartSize;
float *mergeTemp;
using namespace std;
float *readFile(MPI_File *in, const int rank, const int size, const int inputSize){
    if(rank < size){
    int remain = inputSize % size;
    int i;
    partSize = inputSize / size ;
    rightPartSize = partSize + (remain>rank+1);
    startPos = rank * partSize + ((remain>rank)?rank:remain);
    partSize += (remain>rank);

    float *localData = (float *) malloc(sizeof(float) * partSize);

    MPI_File_read_at(*in,sizeof(float) * startPos, localData, partSize, MPI_FLOAT, MPI_STATUS_IGNORE);

    return localData;
    }
    return NULL;
}
int merge(float *localData, float *neiborData, float *mergeTemp, int recvNum){

    int i = 0, j = 0, k = 0;
    int start;
    start = i = k = upper_bound(localData, localData+partSize, neiborData[0]) - localData;

    while(i < partSize && j < recvNum){
        if(neiborData[j] < localData[i]){
            mergeTemp[k] = neiborData[j];
            j++;
        }
        else{
            mergeTemp[k] = localData[i];
            i++;
        } 
        k++;
    }
    memcpy(mergeTemp+k, localData+i, sizeof(float)*(partSize-i));
    memcpy(localData+start, mergeTemp+start, sizeof(float)*(partSize-start));
    
    return 0;
}
void odd_even_sort(float *localData, const int rank, const int size,const int inputSize){
    float *neiborData;
    MPI_Request request1, request2;
    MPI_Status status;
    if(rank < size){
        mergeTemp = (float *)malloc(sizeof(float) * (rightPartSize + partSize));
        neiborData = mergeTemp + partSize;
        sort(localData, localData+partSize);
        
    }
    int running = 1;
    int swapped;
    int recvnum = 0;
    float tmp;
    int i;
    while(running){
        swapped = 0;
        running = 0;

    //even phase:
    if(rank < size){
        if(rank & 1){
            //even phase
            
            MPI_Irecv(&tmp, 1, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request1);
            MPI_Wait(&request1, &status);
            recvnum = lower_bound(localData, localData+partSize, tmp) - localData;
            MPI_Isend(&recvnum, 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &request2);
            MPI_Wait(&request2, &status);
            if(recvnum != 0){
                MPI_Isend(localData, recvnum, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request1);
                MPI_Wait(&request1, &status);
                MPI_Irecv(localData, recvnum, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request2);
                swapped |= 1;
                MPI_Wait(&request2, &status);
                
            }
            
            //odd phase
            if(rank != size - 1){
                tmp = localData[partSize-1];
                MPI_Isend(&tmp, 1, MPI_FLOAT, rank+1, 3,  MPI_COMM_WORLD, &request1);
                MPI_Wait(&request1, &status);
                MPI_Irecv(&recvnum, 1, MPI_INT, rank+1, 3, MPI_COMM_WORLD, &request2);
                MPI_Wait(&request2, &status);
                if(recvnum != 0){

                    MPI_Irecv(neiborData, recvnum, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request1);
                    MPI_Wait(&request1, &status);
                    
                    merge(localData, neiborData, mergeTemp, recvnum);
                    swapped |= recvnum;
                    MPI_Isend(neiborData, recvnum, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request2);
                    MPI_Wait(&request2, &status);

                }
                
            }
        }
        else {
            if( rank != size - 1){
                tmp = localData[partSize-1];
                MPI_Isend(&tmp, 1, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request1);
                MPI_Wait(&request1, &status);
                MPI_Irecv(&recvnum, 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &request2);
                MPI_Wait(&request2, &status);
                if(recvnum != 0){
                    MPI_Irecv(neiborData, recvnum, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request1);
                    MPI_Wait(&request1, &status);
                    merge(localData, neiborData, mergeTemp, recvnum);
                    MPI_Isend(neiborData, recvnum, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request2);
                    MPI_Wait(&request2, &status);
                    swapped |= recvnum;
                }
                    
            }
            if(rank != 0){
                MPI_Irecv(&tmp, 1, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &request1);
                MPI_Wait(&request1, &status);
                recvnum = lower_bound(localData, localData+partSize, tmp) - localData;
                MPI_Isend(&recvnum, 1, MPI_INT, rank-1, 3, MPI_COMM_WORLD, &request2);
                MPI_Wait(&request2, &status);
                if(recvnum != 0){
                    MPI_Isend(localData, recvnum, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request1);
                    MPI_Wait(&request1, &status);
                    MPI_Irecv(localData, recvnum, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request2);
                    
                    swapped |= recvnum;
                    MPI_Wait(&request2, &status);
                    
                }
            }
            
        }
    

        }
        MPI_Allreduce(&swapped, &running, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    
    }



}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    int rank, size;
    int inputSize = atoi(argv[1]);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size > inputSize){
        size = inputSize;
    }
    MPI_File in, out;
    float *localData;

    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    localData = readFile(&in, rank, size, inputSize);

    MPI_File_close(&in);

    odd_even_sort(localData, rank, size, inputSize);

    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &out);

    if(rank < size){
        //MPI_File_set_view(out, sizeof(float)*startPos,MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
        //MPI_File_write_all(out, localData, partSize, MPI_FLOAT, MPI_STATUS_IGNORE);

        MPI_File_write_at(out, sizeof(float)*startPos, localData, partSize, MPI_FLOAT, MPI_STATUS_IGNORE);

        free(localData);
        free(mergeTemp);
    
    }
    
    MPI_File_close(&out);
    MPI_Finalize();
    return 0;
}
