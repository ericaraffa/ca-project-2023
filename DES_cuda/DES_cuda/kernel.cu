
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <thread>

#include <fstream>
#include "matrices.h"
#include "numberConversions.h"

#include <vector>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

//#define N_THREAD 16
#define N_BLOCKS 16
#define THREADS_PER_BLOCK 128
#define DES_BLOCK_SIZE 8
#define ROUNDKEY_SIZE 48
#define ROUNDKEY_NUM 16


void my_strncpy_h(char* dest, const char* src, int dim) {

    for (int i = 0; i < dim; i++) {
        dest[i] = src[i];
    }
}

void my_reverse_h(char* roundKeys, int len_key, int num_key) {

    int dim = len_key*num_key;
    
    char* end_ptr = roundKeys + (num_key-1) * len_key; //Last key
    char* begin_ptr = roundKeys; 
    char* tmp = (char*)malloc(len_key);

    /*
    printf("BEGIN\n");
    for (int j = 0; j < len_key; j++) { 
        printf("%c", *(begin_ptr + j)); 
    }
    printf("\n");

    printf("END\n");
    for (int j = 0; j < len_key; j++) { 
        printf("%c", *(end_ptr + j)); 
    }
    printf("\n");
    */
    

    for (int i = 0; i < (num_key)/2; i++) {

        my_strncpy_h(tmp, begin_ptr, len_key);
        my_strncpy_h(begin_ptr, end_ptr, len_key);
        my_strncpy_h(end_ptr, tmp, len_key); //old begin

        // update pointers positions
        begin_ptr = begin_ptr + len_key;
        end_ptr = end_ptr - len_key;
    }

    free(tmp);
}

__device__ void opt_permutation(char *key, int* arr, int n, char *ret)
{
    for (int i = 0; i < n; i++) {
        ret[i] = key[arr[i] - 1]; 
    }
}

__device__ void shift_left(char *k, int shifts, char *ret)
{
    for (int i = 0; i < 28; i++) {
        ret[i] = k[(i + shifts) % 28];
    }     
}

__device__ void opt_xorOperation(char* str1, char* str2, int n, char* xored)
{

    for (int i = 0; i < n; i++) {
        xored[i] = ((str1[i] - '0') ^ (str2[i] - '0')) + '0';
    }
}

__global__ void encryption(char* pt, char* roundKeys, char* ret_cipher, int pt_size, int groups_per_thread)
{

    extern __shared__ char roundKeys_s[];

    // Control
    long index_th = threadIdx.x * groups_per_thread * DES_BLOCK_SIZE;
    long index_block = blockIdx.x * blockDim.x * groups_per_thread * DES_BLOCK_SIZE;
    long index = index_th + index_block;

    if (index_th == 0) {
        my_strncpy(roundKeys_s, roundKeys, ROUNDKEY_SIZE * ROUNDKEY_NUM * sizeof(char));
    }
    __syncthreads();

    if (index < (pt_size - (groups_per_thread * DES_BLOCK_SIZE))) {

        // String to binary 
        char* pt_bin = (char*)malloc(DES_BLOCK_SIZE * 8);

        // For first permutation
        char* pt_perm = (char*)malloc(64);

        // For splitting 64 bits plain text to LPT and RPT of 32 bits each
        char* left = (char*)malloc(32);
        char* right = (char*)malloc(32);

        // For single roundKey
        char* expandedRPT = (char*)malloc(ROUNDKEY_SIZE);
        char* xored = (char*)malloc(ROUNDKEY_SIZE);

        //s_box_i string array for storing the 4 bits outputs 
        char* s_box_i = (char*)malloc(sizeof(int));

        //s_box_final stores the final result from s box. i.e. concat all the result array elements
        char* s_box_final = (char*)malloc(32);

        // For P-box permutation and xor
        char* p_box_perm = (char*)malloc(32);
        char* p_box_xored = (char*)malloc(32);

        // For final result
        char* cipher_perm = (char*)malloc(64);
        char* cipher = (char*)malloc(DES_BLOCK_SIZE);

        for (int i = 0; i < groups_per_thread; i++) {

            // String to binary 
            str2bin(&pt[index + (i * DES_BLOCK_SIZE)], DES_BLOCK_SIZE, pt_bin);

            // Initial permutation Process
            opt_permutation(pt_bin, initialPermutation, 64, pt_perm);

            // Splitting of 64bits plain text to LPT and RPT of 32 bits each
            my_strncpy(left, pt_perm, 32);
            my_strncpy(right, pt_perm + 32, 32);

            for (int j = 0; j < ROUNDKEY_NUM; j++) {
                // Expansion Permutation 
                opt_permutation(right, dBox, ROUNDKEY_SIZE, expandedRPT);

                // XOR of RoundKey and expandedRPT
                opt_xorOperation(roundKeys_s + (j * ROUNDKEY_SIZE), expandedRPT, ROUNDKEY_SIZE, xored);

                // S-boxes
                for (int k = 0; k < 8; k++) {
                    //the value of '0' is 48, '1' is 49 and so on. but since we are referring the matrix index, we are interested in 0,1,..
                    //So, the '0' should be subtracted . i.e. the 49 value of '1' will be 49-48=1.
                    int row = 2 * int(xored[k * 6]) + int(xored[k * 6 + 5]);
                    int col = 8 * int(xored[k * 6 + 1]) + 4 * int(xored[k * 6 + 2]) + 2 * int(xored[k * 6 + 3]) + int(xored[k * 6 + 4]);

                    int val = sbox[k][row][col];

                    decimalToBinary(val, s_box_i);

                    my_strncpy(s_box_final + (k * sizeof(int)), s_box_i, sizeof(int));
                }

                // P-Box Permutation 
                opt_permutation(s_box_final, pbox, 32, p_box_perm);
                // XOR of left and p_box_perm 
                opt_xorOperation(p_box_perm, left, 32, p_box_xored);
                // Update left
                my_strncpy(left, p_box_xored, 32);
                // Swap left and right in every rounds except the last round
                if (j != 15) {
                    char* tmp;
                    tmp = left;
                    left = right;
                    right = tmp;
                }
            }
            // Left and Right combined
            my_strncpy(pt_perm, left, 32);
            my_strncpy(pt_perm + 32, right, 32);

            // Final Permutation to obtain 64bits cipher text
            opt_permutation(pt_perm, finalPermutation, 64, cipher_perm);

            bin2str(cipher_perm, DES_BLOCK_SIZE * 8, cipher);      

            my_strncpy(&ret_cipher[index + (i * DES_BLOCK_SIZE)], cipher, DES_BLOCK_SIZE);

        }
        //FREE
        free(pt_bin);
        free(pt_perm);
        free(left);
        free(right);
        free(expandedRPT);
        free(xored);
        free(s_box_i);
        free(s_box_final);
        free(p_box_perm);
        free(p_box_xored);
        free(cipher_perm);
        free(cipher);
    }   
}

__global__ void generate_roundKeys(char *key, int key_size, char *ret_roundKeys, int roundKey_num) {

    char * key_bin = (char*)malloc(key_size * 4);
    char * key_perm_shift = (char*)malloc(56);

    convertToBinary(key, key_size, key_bin); 
    opt_permutation(key_bin, keyTransformation, 56, key_perm_shift); // key without parity 

    free(key_bin);
   
    //Splitting 56 bit keys to left and right of 28 bits each
    char *left = (char*)malloc(28);
    char *right = (char*)malloc(28);

    for (int i = 0; i < roundKey_num; i++) {
        //Left Shift and Right Shift done to the respective left and right keys in each round
        shift_left(key_perm_shift, shiftsMatrix[i], left);
        shift_left((key_perm_shift + 28), shiftsMatrix[i], right);
        my_strncpy(key_perm_shift, left, 28);
        my_strncpy(key_perm_shift + 28, right, 28);

        // Key Compression : Converting 56 bit key to 48 bit combined key
        opt_permutation(key_perm_shift, keyCompression, 48, (ret_roundKeys + (i * 48)));
    }

    free(key_perm_shift);
    free(left);
    free(right);
}


int main()
{

    cudaEvent_t cuda_start, cuda_stop;  
    cudaEventCreate(&cuda_start); 
    cudaEventCreate(&cuda_stop); 
    float milliseconds_enc = 0; 
    float milliseconds_dec = 0;

    FILE* file = fopen("pt_10000kb.txt", "rb");
    if (!file) {
        cout << "dead" << endl;
    }

    char *pt;
    long FILE_SIZE = 0;
    long PT_SIZE = 0;
    long PADDING = 0;
    int c, n_groups; 

    fseek(file, 0L, SEEK_END);
    FILE_SIZE = ftell(file);
    rewind(file);

    while ((FILE_SIZE + PADDING) % DES_BLOCK_SIZE != 0) {
        PADDING++; 
    }

    n_groups = (FILE_SIZE + PADDING) / DES_BLOCK_SIZE;

    //Allocate mem for pt
    PT_SIZE = n_groups * DES_BLOCK_SIZE * sizeof(char);
    pt = (char*)malloc(PT_SIZE);

    //Copy the file
    for (int i = 0; i < n_groups; i++) {
        if (i == n_groups - 1) { 
            fread(pt + (i * DES_BLOCK_SIZE), DES_BLOCK_SIZE - PADDING, 1, file);

            for (int j = (i * DES_BLOCK_SIZE) + DES_BLOCK_SIZE - PADDING; j < n_groups * DES_BLOCK_SIZE; j++) {
                pt[j] = ' ';
            }
        }
        else {
            fread(pt + (i * DES_BLOCK_SIZE), DES_BLOCK_SIZE, 1, file);
        }
    }
    fclose(file);
    
    // GENERATE KEYS
    const int key_size = 17;
    char key[key_size] = "ABC12532110EDA56";
    char* roundKeys;

    roundKeys = (char*)malloc(ROUNDKEY_SIZE * ROUNDKEY_NUM * sizeof(char));


    // GPU
    char* d_roundKeys;
    char* d_key;

    cudaMalloc((void**)&d_key, key_size); 
    cudaMalloc((void**)&d_roundKeys, ROUNDKEY_SIZE * ROUNDKEY_NUM * sizeof(char));

    cudaMemcpy(d_key, &key, key_size, cudaMemcpyHostToDevice);

    generate_roundKeys<<<1,1>>>(d_key, key_size, d_roundKeys, ROUNDKEY_NUM);
    cudaDeviceSynchronize();

    cudaFree(d_key);

    cudaMemcpy(roundKeys, d_roundKeys, ROUNDKEY_SIZE * ROUNDKEY_NUM *sizeof(char), cudaMemcpyDeviceToHost);

    //	ENCRYPTION: creation of the vector cipher (of size pt.size()) and a vector of threads. After defining the N_THREAD the 'blocks'
    //				of the plain text are splitted by the num of threads.
    //				The function thread t() create a thread and it will be pushed inside the array threads.


    char* cipher = (char*)malloc(PT_SIZE);

    // GPU
    char* d_pt;
    char* d_cipher;

    // Allocate
    cudaMalloc((void**)&d_pt, PT_SIZE);
    cudaMalloc((void**)&d_cipher, PT_SIZE);
    // Copy inputs
    cudaMemcpy(d_pt, pt, PT_SIZE, cudaMemcpyHostToDevice);

    dim3 grid(N_BLOCKS);
    dim3 block(THREADS_PER_BLOCK);

    int groups_per_thread = n_groups / (block.x * grid.x);

    int shared_mem_dim = ROUNDKEY_SIZE * ROUNDKEY_NUM * sizeof(char);

    cudaEventRecord(cuda_start); 
    encryption <<<grid, block, shared_mem_dim>>> (d_pt, d_roundKeys, d_cipher, PT_SIZE, groups_per_thread);
    cudaEventRecord(cuda_stop);
 
    // Copy results
    cudaMemcpy(cipher, d_cipher, PT_SIZE, cudaMemcpyDeviceToHost);

    
    cudaEventSynchronize(cuda_stop); 
    cudaEventElapsedTime(&milliseconds_enc, cuda_start, cuda_stop); 
    
    // Decryption : Reversing the round keys and executing the encryption process to get Plain Text

    my_reverse_h(roundKeys, ROUNDKEY_SIZE, ROUNDKEY_NUM); //ok

    // GPU decrypt
    // Copy input
    cudaMemcpy(d_roundKeys, roundKeys, ROUNDKEY_SIZE * ROUNDKEY_NUM * sizeof(char), cudaMemcpyHostToDevice); 
    cudaMemcpy(d_cipher, cipher, PT_SIZE, cudaMemcpyHostToDevice); 

    cudaEventRecord(cuda_start);
    encryption <<<grid, block, shared_mem_dim>>> (d_cipher, d_roundKeys, d_pt, PT_SIZE, groups_per_thread);
    cudaEventRecord(cuda_stop);

    // Copy results
    cudaMemcpy(pt, d_pt, PT_SIZE, cudaMemcpyDeviceToHost); 
    
    cudaEventSynchronize(cuda_stop);
    cudaEventElapsedTime(&milliseconds_dec, cuda_start, cuda_stop);   

    // Cleanup
    cudaFree(d_pt); 
    cudaFree(d_roundKeys); 
    cudaFree(d_cipher);
    
    free(pt);
    free(roundKeys);
    free(cipher);

    cout << "Time taken (encryption+decryption): " << milliseconds_enc + milliseconds_dec << " ms" << endl;

}
