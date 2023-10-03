
#include <omp.h>
#include <thread>

#include <fstream>
#include "matrices.h"
#include "numberConversions.h"

#include <vector>
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

#define N_THREAD 16

string permutation(string str, int* arr, int n)
{
	string res = "";
	for (int i = 0; i < n; i++) {
		res += str[arr[i] - 1];
	}
	return res;
}

string opt_permutation(string str, int* arr, int n)
{
	string result;
	result.reserve(n);

	for (int i = 0; i < n; i++) {
		result += str[arr[i] - 1];
	}
	return result;
}

string shift_left(string k, int shifts)
{
	string s = "";
	for (int i = 0; i < shifts; i++) {
		for (int j = 1; j < 28; j++) {
			s += k[j];
		}
		s += k[0];
		k = s;
		s = "";
	}
	return k;
}

string xorOperation(string str1, string str2)
{
	string xored = "";
	for (int i = 0; i < str1.size(); i++) {
		if (str1[i] == str2[i])
			xored += "0";
		else
			xored += "1";
	}
	return xored;
}

string opt_xorOperation(string str1, string str2)
{
	string xored;
	xored.reserve(str1.size());

	for (int i = 0; i < str1.size(); i++) {
		xored += ((str1[i] - '0') ^ (str2[i] - '0')) + '0';
	}
	return xored;
}

string encryption(string pt, vector<string> roundKeys)
{
	// String to binary 
	pt = str2bin(pt);
	// Initial permutation Process
	pt = opt_permutation(pt, initialPermutation, 64);
	// Splitting of 64bits plain text to LPT and RPT of 32 bits each
	string left = pt.substr(0, 32);
	string right = pt.substr(32, 32);
	for (int i = 0; i < 16; i++) {
		// Exapansion Permutation 
		string expandedRPT = opt_permutation(right, dBox, 48);
		// XOR of RoundKey and expandedRPT
		string x = opt_xorOperation(roundKeys[i], expandedRPT);
		// S-boxes
		//result string array for storing the 4 bits outputs for 8*6 bits input
		string result[8];
		//res stores the final result from s box. i.e. concat all the result array elements
		string res = "";
		for (int i = 0; i < 8; i++) {
			//the value of '0' is 48, '1' is 49 and so on. but since we are referring the matrix index, we are interested in 0,1,..
			//So, the '0' should be subtracted . i.e. the 49 value of '1' will be 49-48=1.
			int row = 2 * int(x[i * 6] - '0') + int(x[i * 6 + 5] - '0');
			int col = 8 * int(x[i * 6 + 1] - '0') + 4 * int(x[i * 6 + 2] - '0') + 2 * int(x[i * 6 + 3] - '0') + int(x[i * 6 + 4] - '0');
			int val = sbox[i][row][col];
			result[i] = decimalToBinary(val);
		}
		for (int i = 0; i < 8; i++)
			res += result[i];
		// P-Box Permutation 
		res = opt_permutation(res, pbox, 32);
		// XOR of left and res 
		x = opt_xorOperation(res, left);
		left = x;
		// Swap left and right in every rounds except the last round
		if (i != 15) {
			swap(left, right);
		}
	}
	// Left and Right combined
	string combined = left + right;
	// Final Permutation to obtain 64bits cipher text
	string cipher = bin2str(opt_permutation(combined, finalPermutation, 64));
	return cipher;
}

// void encryption_split(int start, int end, vector<string>& cipher, vector<string> pt, vector<string> roundKeys, int & throughput) {
void encryption_split(int start, int end, vector<string>& cipher, vector<string> pt, vector<string> roundKeys) {

	for (int i = start; i < end; i++) {
		cipher[i] = encryption(pt[i], roundKeys);

		// THROUGHPUT
		//throughput += 8;
	}

}


int main() {

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	//Reading from the file and storing the plaintext in pt vector
	vector<string> pt;
	string c, temp;
	ifstream MyReadFile("pt_10000kb.txt"); //9kb
	while (getline(MyReadFile, c)) {
		temp += c;
	}
	while (temp.length() % 8 != 0) {
		temp += " ";
	}
	for (int i = 0; i < temp.length(); i = i + 8) {
		pt.emplace_back(temp.substr(i, 8));
	}
	MyReadFile.close();
	string key = "ABC12532110EDA56";
	key = convertToBinary(key);
	key = opt_permutation(key, keyTransformation, 56); // key without parity 
	//Splitting 56 bit keys to left and right of 28 bits each
	string left = key.substr(0, 28);
	string right = key.substr(28, 28);
	vector<string> roundKeys(16); // Declaring vector for storing keys of 16 rounds
	for (int i = 0; i < 16; i++) {
		//Left Shift and Right Shift done to the respective left and right keys in each round
		left = shift_left(left, shiftsMatrix[i]);
		right = shift_left(right, shiftsMatrix[i]);
		string combinedkey = left + right;
		// Key Compression : Converting 56 bit key to 48 bit combined key
		string RoundKey = opt_permutation(combinedkey, keyCompression, 48);
		roundKeys.emplace_back(RoundKey);
	}


	//	ENCRYPTION: creation of the vector cipher (of size pt.size()) and a vector of threads. After defining the N_THREAD the 'blocks'
	//				of the plain text are splitted by the num of threads.
	//				The function thread t() create a thread and it will be pushed inside the array threads.


	vector<string> ciphertext;
	vector<string> cipher(pt.size());
	vector<thread> threads;

	//vector<int> throughput(N_THREAD, 0);

	int split_th = pt.size() / N_THREAD;

	std::chrono::steady_clock::time_point begin_encryption = std::chrono::steady_clock::now();

	for (int j = 0; j < N_THREAD; j++) {
		int start = j * split_th;
		int end = j == N_THREAD - 1 ? pt.size() : split_th * (j + 1); // if j is last thread, takes the leftovers 

		//THROUGHPUT
		//thread t(encryption_split, start, end, ref(cipher), pt, roundKeys, ref(throughput[j]));
		thread t(encryption_split, start, end, ref(cipher), pt, roundKeys);

		threads.emplace_back(move(t));
	}

	// THROUGHPUT 200ms 
	//std::this_thread::sleep_for(std::chrono::milliseconds(200)); 

	// Now wait for the threads to finish,
	// We need to wait otherwise main thread might reach an end before the multiple threads finish their work
	for (auto& t : threads) {
		t.join();
		
		// THROUGHPUT
		//t.detach();
		//t.~thread();
	}
	std::chrono::steady_clock::time_point end_encryption = std::chrono::steady_clock::now();

	/*
	ofstream throughput_out("throughput.txt", ios::app); 
	throughput_out << "#thread: " << N_THREAD << endl;
	for (int i = 0; i < N_THREAD; i++) {
		throughput_out << throughput[i] << " ";
	}
	throughput_out << endl;

	throughput_out.close();
	*/

	//cout << "cipher:" << cipher.size() << endl;


	for (int i = 0; i < cipher.size(); i++) {
		ciphertext.emplace_back(cipher[i]);
	}

	//Writing Cipher Text to File
	ofstream writeObj;
	remove("encrypted.txt");
	writeObj.open("encrypted.txt", ios::app);

	if (!writeObj)
	{
		return 0;
	}


	for (int j = 0; j < ciphertext.size(); j++) {
		writeObj << ciphertext[j];
	}
	writeObj.close();

	//Decryption : Reversing the round keys and executing the encryption process to get Plain Text
	reverse(roundKeys.begin(), roundKeys.end());
	string decrypted;
	vector <string> text(ciphertext.size());

	split_th = ciphertext.size() / N_THREAD;

	threads.clear();

	std::chrono::steady_clock::time_point begin_decryption = std::chrono::steady_clock::now();

	for (int j = 0; j < N_THREAD; j++) {
		int start = j * split_th;
		int end = j == N_THREAD - 1 ? ciphertext.size() : split_th * (j + 1); // if j is last thread, takes the leftovers 

		thread t(encryption_split, start, end, ref(text), ciphertext, roundKeys);

		threads.emplace_back(move(t));
	}

	for (auto& t : threads) {
		t.join();
	}

	std::chrono::steady_clock::time_point end_decryption = std::chrono::steady_clock::now();

	for (int i = 0; i < text.size(); i++) {
		decrypted += text[i];
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	cout << "Time taken (main): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " ms" << endl;
	cout << "Time taken (encryption+decryption): " << std::chrono::duration_cast<std::chrono::milliseconds>(end_encryption - begin_encryption +
		end_decryption - begin_decryption).count() << " ms" << endl;
	
}