#include<mpi.h>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include <algorithm>
#include<ctime>
#include<stack>
#define LEFT  0
#define RIGHT  1
#define MIN 0
#define MAX 1

long int localSize;
long int partnerSize[2];
int partnerRank[2];
float *mergeAry[2];
float *localData;
float *partnerData;
int cur;
int startPos;

using namespace std;

void Init(MPI_File *in, const int rank, const int size, const int inputSize){
    if(rank < size){
    int remain = inputSize % size;
    localSize = inputSize / size ;

    partnerRank[LEFT] = rank-1;
    partnerRank[RIGHT] = rank+1;
    partnerSize[LEFT] = localSize + (remain > partnerRank[LEFT]);
    partnerSize[RIGHT] = localSize + (remain > partnerRank[RIGHT]);
    startPos = rank * localSize + ((remain>rank)?rank:remain);
    localSize += (remain>rank);
    cur = 0;
    mergeAry[0] = new float[localSize+partnerSize[LEFT]];
    mergeAry[1] = new float[localSize+partnerSize[LEFT]];
    localData = mergeAry[cur];
    partnerData = mergeAry[cur]+localSize;
    //MPI_File_read_at(*in,sizeof(float) * startPos, localData, localSize, MPI_FLOAT, MPI_STATUS_IGNORE); 
    MPI_File_read_shared(*in, localData, localSize, MPI_FLOAT, MPI_STATUS_IGNORE);
    sort(localData, localData + localSize);
    }
    
}
int Exchange(int rank, int size, int dir, int mode){
    if(rank < 0 || rank >= size || partnerRank[dir] < 0 || partnerRank[dir] >= size) return 0;
    int sendnum, recvnum, offset = 0;
    float tmp, ntmp;
    tmp = (mode == MIN) ? localData[localSize-1] : localData[0];
    MPI_Sendrecv(&tmp, 1, MPI_FLOAT, partnerRank[dir], mode,
                        &ntmp, 1, MPI_FLOAT, partnerRank[dir], 1-mode, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(mode == MIN){
        sendnum = localSize - (upper_bound(localData, localData+localSize, ntmp) - localData);
        offset = localSize - sendnum;    
    }
    else
        sendnum = lower_bound(localData, localData+localSize, ntmp) - localData;
    
    MPI_Sendrecv(&sendnum, 1, MPI_INT, partnerRank[dir], mode,
                &recvnum, 1, MPI_INT, partnerRank[dir], 1-mode, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    if(recvnum == 0) return 0;
    MPI_Sendrecv(localData+offset, sendnum, MPI_FLOAT, partnerRank[dir], mode,
                 partnerData, recvnum, MPI_FLOAT, partnerRank[dir], 1-mode,
                  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int next = 1-cur;
    if(mode == MIN){
        merge(localData+offset, localData+localSize,
            partnerData, partnerData+recvnum,
            mergeAry[next]+offset);
        memcpy(mergeAry[next], localData, sizeof(float)*offset);
        localData = mergeAry[next];
        partnerData = mergeAry[next] + localSize;
    }
    else{
        merge(localData, localData+sendnum,
            partnerData, partnerData+recvnum,
            mergeAry[next]+partnerSize[dir]-recvnum);
        memcpy(mergeAry[next]+partnerSize[dir]+sendnum, localData+sendnum, sizeof(float)*(localSize-sendnum));
        localData = mergeAry[next] +partnerSize[dir];
        partnerData = mergeAry[next];
    }

    cur = next;
    return 1;
}
void odd_even_sort(int rank, int size){
    
        int isSwap , running=1;
        while(running){
            //even phase
            isSwap = 0;
            if(rank & 1)
                isSwap |= Exchange(rank, size, LEFT, MAX);
            else
                isSwap |= Exchange(rank, size, RIGHT, MIN);
            
            //odd phase
            if(rank & 1)
                isSwap |= Exchange(rank, size, RIGHT, MIN);
            else
                isSwap |= Exchange(rank, size, LEFT, MAX);
            MPI_Allreduce(&isSwap, &running, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
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

    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    Init(&in, rank, size, inputSize);
    MPI_File_close(&in);

    odd_even_sort(rank, size);
    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &out);
    if(rank < size){
        
        MPI_File_write_at(out, sizeof(float)*startPos, localData, localSize, MPI_FLOAT, MPI_STATUS_IGNORE);
        
        delete [] mergeAry[0];
        delete [] mergeAry[1];
    }
    MPI_File_close(&out);
    
    MPI_Finalize();
    return 0;
}
