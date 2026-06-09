/// @file helper.hpp

// https://groups.google.com/g/tmb-users/c/xg3pQwNbxEY/m/B0lGUWG3AQAJ
namespace CppAD {
  void PrintFor(const char* before, const double& var) {  }
}

// List of Type vectors
template<class Type>
struct LOV_t : vector<vector<Type> > {
  LOV_t(SEXP x){ /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP v = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(v);
    }
  }
};

// List of int vectors
template<class Type>
struct LOVI_t : vector<vector<int> > {
  LOVI_t(SEXP x){ /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP v = VECTOR_ELT(x, i);
      (*this)(i) = asVector<int>(v);
    }
  }
};

// List of dense matrices
template<class Type>
struct LODM_t : vector<matrix<Type> > {
  LODM_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP dm = VECTOR_ELT(x, i);
      if(!isMatrix(dm))
        error("Not a dense matrix");
      (*this)(i) = asMatrix<Type>(dm);
    }
  }
};

// List of sparse matrices (https://github.com/kaskr/adcomp/issues/96#issuecomment-109792130)
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!isValidSparseMatrix(sm))
        error("Not a sparse matrix");
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};

template<class Type>
void fillSparseMatId(SparseMatrix<Type>& Id)
{
  Id.setZero();
  for(int i=0; i< Id.rows(); ++i)
    Id.coeffRef(i,i) = 1.0;
}

template<class Type>
void invvec(const vector<Type>& x,
            const int& ncol,
            matrix<Type>& mat)
{
  int n = x.size();
  int nrow = (int) n / ncol;
  int idx = 0;
  for (int idx_col = 0; idx_col < ncol; ++idx_col) {
    mat.col(idx_col) = x.segment(idx, nrow);
    idx = idx + nrow;
  }
}

template<class Type>
void invvecMixes(const vector<Type>& x,
                 const int& nb_comps,
                 matrix<Type>& mat)
{
  int n = x.size();
  int nrow = (int) n / nb_comps;
  for (int j = 0; j < nb_comps; ++j) {
    for (int i = 0; i < nrow; ++i) {
      mat(i, j) = x[i * nb_comps + j];
    }
  }
}

template<class Type>
void invvecMixes(const vector<Type>& x,
                 const int& nb_comps,
                 array<Type>& arr)
{
  int n = x.size();
  int nrow = (int) n / nb_comps;
  for (int j = 0; j < nb_comps; ++j) {
    for (int i = 0; i < nrow; ++i) {
      arr(i, j) = x[i * nb_comps + j];
    }
  }
}
