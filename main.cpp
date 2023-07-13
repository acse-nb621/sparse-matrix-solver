#include <iostream> /* cout */
#include "sparseMatrix.h"
#include "Matrix.h"
#include <algorithm> 
#include <random>
#include <time.h> /* time */
#include <stdio.h> /*NULL*/
#include <stdlib.h> /* srand, rand */
using namespace std;


template <class A>
//templated wrapper function to solve a linear system by either Jacobi ot the Gauss-Seidel method
vector<double> solve(sparseMatrix<A>& mat, vector<double>& vec, char s = 'j');

//generate a random tri-diagonal, symmetric, diagonally dominant matrix of doubles
sparseMatrix<int> generate_sparse_int(int size);

//generate a random vector of double
vector<double> generate_vector(int size);

template <class A>
//templated function to print a vector
void printv(vector<A>& vec);



int main(){
    vector<int> sizes = {10,50,100};
    for(int& i: sizes){
        cout<<"\nBelow is a randomly generated "<<i<<"x"<<i<<" sparse matrix of doubles.\nThis matrix is symmetric, tri-diagonal and diagonally dominant.\n";
        sparseMatrix<int> mat = generate_sparse_int(i);
        mat.printMatrix();
        cout<<"\nBelow is a randomly generated "<<i<<"x"<<i<<" vector of doubles.\n";
        vector<double> vec = generate_vector(i);
        cout<<endl;
        printv(vec);
        cout<<"\nBelow is the solution obtained using the Jacobi Method.\n";
        vector<double> solsj = solve(mat,vec,'j');
        cout<<endl;
        printv(solsj);
        cout<<endl;
        cout<<"\nBelow is the solution obtained using the Gauss-Seidel Method.\n";
        vector<double> solsg = solve(mat,vec,'g');
        cout<<endl;
        printv(solsg);
        cout<<endl;
    }
}

sparseMatrix<int> generate_sparse_int(int size){
    srand(time(NULL));
    int* mat = new int[size*size]();
    for(int i = 0; i < size*(size-1) ; i += size){
        
        mat[i + i/size] = (int)rand()/(RAND_MAX)+1 +(rand()%(5-1)) + 5;
        mat[i + i/size + 1] = (int)rand()/(RAND_MAX)+1 +(rand()%(5-1));
        mat[i + i/size + size] = mat[i + i/size + 1];
    }
    // fix last elemt
    mat[size*size-1] = (int)rand()/(RAND_MAX)+1 +(rand()%(5-1)) + 5;
    sparseMatrix<int> spm(size,size,mat);
    return spm;
}

template <class A>
void printv(vector<A>& vec){
    for(int i = 0; i < vec.size(); i++){
        cout<<vec[i]<<" ";
    }
    cout<<endl;
}

template <class B>
vector<double> solve(sparseMatrix<B>& mat, vector<double>& vec, char s){
    vector<double> ans;
    if(s=='j'){ans = mat.solve_jacobi(vec);}
    else if(s=='g'){ans = mat.solve_gauss(vec);}
    return ans;
}

vector<double> generate_vector(int size){
    srand(time(NULL));
    vector<double> vec;
    for(int i = 0; i<size; i++){
        vec.push_back((double)rand()/(RAND_MAX)+1 +(rand()%9));
    }
    return vec;
}