#pragma once
#include<fstream>
template <class T>
class Matrix
{
public:

   // constructor where we want to preallocate ourselves
   Matrix(int rows, int cols, bool preallocate);
   // constructor where we already have allocated memory outside
   Matrix(int rows, int cols, T *values_ptr);
   // destructor
   virtual ~Matrix();
   // Print out the values in our matrix
   void printValues();
   // Print out the matrix
   virtual void printMatrix();
   void printMatrixToFile(std::ofstream& out); 
   // Explicitly using the C++11 nullptr here
   T *values = nullptr;   
   int rows = -1;
   int cols = -1;
   
protected:
   int size_of_values = -1;
   bool preallocated = false;
   
};
