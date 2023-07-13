#include <vector>
#include "Matrix.h"
#include <cmath>
using namespace std;

//templated class reresenting a sparse matrix in CSR format
template <class T>
class sparseMatrix: public Matrix<T>{
    private:
        //a flag to know whether or ot the values where set by the user
        bool filled = false;

        //the non-zero values in the matrix
        vector<T> sparse_vals;

        //the column of every stored none-zero value
        vector<int> cols_ind; 

        //where every row begins in cols_ind & sparse_vals
        vector<int> row_ptr;

        //private method used to find the row of the element given its index in sparse_vals/cols_ind
        int find_row(int ind);
    public:
        //constructor for a sparse matrix given the row, columns, and values
        sparseMatrix(int rows_, int cols_, T* vals_);

        //constructor for an empty sparse matrix
        sparseMatrix(int rows, int cols);

        //return size of matrix
        int sparse_size(){return this->size_of_values-1;}

        //return the vector of row_ptr
        vector<int> get_row_ptr(){ return row_ptr;}

        //return the vector of cols_ind
        vector<int> get_cols_ind(){ return cols_ind;}

        //return the vector of sparse_vals
        vector<T> get_sparse_vals(){ return sparse_vals;}

        //set the values of an empty sparseMatrix
        void set_matrix_vals(T* arr);
        
        //return the solution of solving a system using the Jacobi method
        vector<double> solve_jacobi(vector<double>& vec);

        //return the solution of solving a system using the Gauss-Seidel method
        vector<double> solve_gauss(vector<double>& vec);

        //return the result of multiplying a sparseMatrix with a vector
        vector<double> mult_sparse_vec(vector<double>& vec);   
};