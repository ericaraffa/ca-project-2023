
#include <string.h>
#include <bitset>
#define numberConversions
using namespace std;

__device__ void my_strncpy(char* dest, const char* src, int dim) {

	for (int i = 0; i < dim; i++) {
		dest[i] = src[i];
	}
}

__device__ void my_swap(char* src, char* dest) {
	char* tmp;
	tmp = src; 
	src = dest; 
	dest = tmp;

}

__device__ void my_strncat(char* src, char* dest, int dim) {

	for (int i = 0; i < dim; i++) {
		dest += src[i];
	}
}

__device__ void my_substr(char* src, int start, int end, char* dest) {

	for (int i = start; i < end; i++) {
		dest[i] = src[i];
	}

}

__device__ void str2bin(char* s, int size, char* dest)
{
	char* str;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 8; j++) {
			dest[j + (i*8)] = (s[i] >> (7 - j)) & 1;
		}		
	}
}

__device__ void binaryToDecimal(char* val, int* res) {

	for (int i = 0; i < 8; i++) {
		*res = *res + int(val[7 - i]) * (1 << i);
	}

}

__device__ void bin2str(char* str, int size, char* dest)
{
	char bin2int[8];
	int res = 0;

	for (int i = 0; i < size; i = i + 8) {
		my_strncpy(bin2int, str+i, 8);
		binaryToDecimal(bin2int, &res);

		char c = (char)res;
		dest[i/8] = c;

		res = 0;
	}
}

__device__ void decimalToBinary(int val, char* res) {
	//When ASCII encoding is used, the integer value of '0' is 48. So, adding '0' i.e. 48 to get the integer value of 1 as '1'.
	//i.e. Converting digit to its corresponding character.
	res[0] = char(val / 8 + '0');
	val = val % 8;
	res[1] = char(val / 4 + '0');
	val = val % 4;
	res[2] = char(val / 2 + '0');
	val = val % 2;
	res[3] = char(val + '0');
}

__device__ void convertToBinary(char* s, int dim, char* res)
{
	for (int i = 0; i < dim - 1; i++) {
		switch (s[i]) {
		case '0':
			my_strncpy(res + (i * 4), "0000", 4);
			break;
		case '1':
			my_strncpy(res + (i * 4), "0001", 4);
			break;
		case '2':
			my_strncpy(res + (i * 4), "0010", 4);
			break;
		case '3':
			my_strncpy(res + (i * 4), "0011", 4);
			break;
		case '4':
			my_strncpy(res + (i * 4), "0100", 4);
			break;
		case '5':
			my_strncpy(res + (i * 4), "0101", 4);
			break;
		case '6':
			my_strncpy(res + (i * 4), "0110", 4);
			break;
		case '7':
			my_strncpy(res + (i * 4), "0111", 4);
			break;
		case '8':
			my_strncpy(res + (i * 4), "1000", 4);
			break;
		case '9':
			my_strncpy(res + (i * 4), "1001", 4);
			break;
		case 'A':
		case 'a':
			my_strncpy(res + (i * 4), "1010", 4);
			break;
		case 'B':
		case 'b':
			my_strncpy(res + (i * 4), "1011", 4);
			break;
		case 'C':
		case 'c':
			my_strncpy(res + (i * 4), "1100", 4);
			break;
		case 'D':
		case 'd':
			my_strncpy(res + (i * 4), "1101", 4);
			break;
		case 'E':
		case 'e':
			my_strncpy(res + (i * 4), "1110", 4);
			break;
		case 'F':
		case 'f':
			my_strncpy(res + (i * 4), "1111", 4);
			break;
		}
	}
}