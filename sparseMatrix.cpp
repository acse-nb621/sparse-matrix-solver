#include <vector>
#include <iostream>
#include <algorithm>
#include "sparseMatrix.h"
#include <string>
using namespace std;

//Error class
class Error{}; 

//helper function to throw error when user tries to operate on empty sparseMatrix
void throw_empty_error(){
    cout<<"This matrix has not been filled with values.\nUse the set_matrix_vals function to fill it.\n";
    throw Error();
}

//helper function to throw error when user tries to fill an already filled sparseMatrix
void throw_filled_error(){
    cout<<"This matrix has already been filled with values.\n";
    throw Error();
}

template <class T>
sparseMatrix<T>::sparseMatrix(int rows, int cols):Matrix<T>(rows,cols,true){}

template <class T>
sparseMatrix<T>::sparseMatrix(int rows_, int cols_, T* vals_):Matrix<T>(rows_,cols_,vals_){
    filled = true;
    row_ptr.push_back(0);
    int row_count{0};
    for(int i = 0; i < this->size_of_values; i++){
        int temp_row = i / this->rows;
        int temp_col = i % this->cols;
        if(this->values[i] != 0){
            row_count++;
            sparse_vals.push_back(this->values[i]);
            cols_ind.push_back(temp_col);
        }
        if(temp_col % this->cols == this->cols - 1){
            row_ptr.push_back(row_ptr.back()+row_count);
            row_count = 0;
        }
    }
}

template<class T>
void sparseMatrix<T>::set_matrix_vals(T* arr){
    if(filled){throw_filled_error();}
    for(int i = 0; i < this->size_of_values; i ++){
        this->values[i] = arr[i];
    }
    filled = true;
    row_ptr.push_back(0);
    int row_count{0};
    for(int i = 0; i < this->size_of_values; i++){
        int temp_row = i / this->rows;
        int temp_col = i % this->cols;
        if(this->values[i] != 0){
            row_count++;
            sparse_vals.push_back(this->values[i]);
            cols_ind.push_back(temp_col);
        }
        if(temp_col % this->cols == this->cols - 1){
            row_ptr.push_back(row_ptr.back()+row_count);
            row_count = 0;
        }
    }

}

template<class T>
vector<double> sparseMatrix<T>::mult_sparse_vec(vector<double>& vec){
    if(!filled){throw_empty_error();}
    vector<double> ans(vec.size(),0);
    for(int i = 0; i < ans.size(); i++){
        for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
            ans[i] += (double)sparse_vals[j]*vec[(double)cols_ind[j]];
        }
    }
    return ans;
}

template <class T>
vector<double> sparseMatrix<T>::solve_gauss(vector<double>& vec){
    if(!filled){throw_empty_error();}
    int size = this->rows;
    vector<double> sols(size,0);
    int itr_num{0};
    while(true){
        itr_num++;
        for(int i = 0; i < size; i++){
            double divider;
            sols[i] = vec[i];
            for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
                if(cols_ind[j] != i){sols[i] -= sparse_vals[j]*sols[cols_ind[j]];}
                else{divider = sparse_vals[j];}
            }
            sols[i] /= divider;
        }
        double err{0};
        vector<double> temp = this->mult_sparse_vec(sols);
        for(int i = 0; i < size; i++){
            err += std::abs(vec[i]-temp[i]);
        }
        if(err <= 0.5 | itr_num > 4999){
            cout<<"Gauss-Seidel method terminated in "<<itr_num<<" iterations, with an error threshold of 0.5, and absolute error of: "<<err<<" .\n";
            break;
        }

    }

    return sols;
}

template<class T>
vector<double> sparseMatrix<T>::solve_jacobi(vector<double>& vec){
    if(!filled){throw_empty_error();}
    int size = this->rows;
    vector<double> sols(size,0);
    int itr_num{0};
    while(true){
        itr_num++;
        vector<double> new_sols(size,0);
        //loop over sols and update
        for(int i = 0; i < size; i++){
            // cout<<"i: "<<i<<endl;
            T divider;
            new_sols[i] = vec[i];
            for(int j = row_ptr[i]; j < row_ptr[i+1]; j++){
                // cout<<"j: "<<j<<endl;
                if(cols_ind[j] != i){new_sols[i] -= sparse_vals[j]*sols[cols_ind[j]];}
                else{divider = sparse_vals[j];}
            }
            new_sols[i] /= divider;
        }
        double err{0};
        vector<double> temp = this->mult_sparse_vec(new_sols);
        for(int i = 0; i < size; i++){
            err += std::abs(vec[i]-temp[i]);
        }
        sols = new_sols;
        if(err <= 0.5 | itr_num > 4999){
            cout<<"Jacobi method terminated in "<<itr_num<<" iterations, with an error threshold of 0.5, and absolute error of: "<<err<<" .\n";
            break;
        }
    }
    return sols;
}


template<class T>
int sparseMatrix<T>::find_row(int ind){
    if(!filled){throw_empty_error();}
    int row{0};
    for(int i = 0; i < this->row_ptr.size() ; i++){
        if(ind+1 >= this->row_ptr[i] && ind+1 <= this->row_ptr[i+1]){
            row = i;
            break;
        }
    }
    return row;
}

template class sparseMatrix<double>;
template class sparseMatrix<int>;
