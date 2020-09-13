//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// p0_starter.h
//
// Identification: src/include/primer/p0_starter.h
//
// Copyright (c) 2015-2020, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <memory>

namespace bustub {

/*
 * The base class defining a Matrix
 */
template <typename T>
class Matrix {
 protected:
  // TODO(P0): Add implementation
  Matrix(int r, int c) {
    rows = r, cols = c;
    linear = new T[rows * cols]();
  }

  // # of rows in the matrix
  int rows;
  // # of Columns in the matrix
  int cols;
  // Flattened array containing the elements of the matrix
  // TODO(P0) : Allocate the array in the constructor. Don't forget to free up
  // the array in the destructor.
  T *linear;

 public:
  // Return the # of rows in the matrix
  virtual int GetRows() = 0;

  // Return the # of columns in the matrix
  virtual int GetColumns() = 0;

  // Return the (i,j)th  matrix element
  virtual T GetElem(int i, int j) = 0;

  // Sets the (i,j)th  matrix element to val
  virtual void SetElem(int i, int j, T val) = 0;

  // Sets the matrix elements based on the array arr
  virtual void MatImport(T *arr) = 0;

  // TODO(P0): Add implementation
  virtual ~Matrix() { delete[] linear; }
};

template <typename T>
class RowMatrix : public Matrix<T> {
 public:
  // TODO(P0): Add implementation
  RowMatrix(int r, int c) : Matrix<T>(r, c) {
    data_ = new T *[r];
    for (int i = 0; i < r; ++i) {
      data_[i] = this->linear + i * c;
    }
  }

  // TODO(P0): Add implementation
  int GetRows() override { return this->rows; }

  // TODO(P0): Add implementation
  int GetColumns() override { return this->cols; }

  // TODO(P0): Add implementation
  T GetElem(int i, int j) override { return data_[i][j]; }

  // TODO(P0): Add implementation
  void SetElem(int i, int j, T val) override { data_[i][j] = val; }

  // TODO(P0): Add implementation
  void MatImport(T *arr) override { memcpy(this->linear, arr, sizeof(T) * this->rows * this->cols); }

  // TODO(P0): Add implementation
  ~RowMatrix() override { delete[] data_; }

 private:
  // 2D array containing the elements of the matrix in row-major format
  // TODO(P0): Allocate the array of row pointers in the constructor. Use these pointers
  // to point to corresponding elements of the 'linear' array.
  // Don't forget to free up the array in the destructor.
  T **data_;
};

template <typename T>
class RowMatrixOperations {
 public:
  // Compute (mat1 + mat2) and return the result.
  // Return nullptr if dimensions mismatch for input matrices.
  static std::unique_ptr<RowMatrix<T>> AddMatrices(std::unique_ptr<RowMatrix<T>> mat1,
                                                   std::unique_ptr<RowMatrix<T>> mat2) {
    // TODO(P0): Add code
    int m1_col = mat1->GetColumns();
    int m1_row = mat1->GetRows();
    int m2_col = mat2->GetColumns();
    int m2_row = mat2->GetRows();
    if (m1_col != m2_col || m1_row != m2_row) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }

    auto mat_ptr = new RowMatrix<T>(m1_col, m1_row);
    for (int i = 0; i < m1_row; ++i) {
      for (int j = 0; j < m1_col; ++j) {
        mat_ptr->SetElem(i, j, mat1->GetElem(i, j) + mat2->GetElem(i, j));
      }
    }
    return std::unique_ptr<RowMatrix<T>>(mat_ptr);
  }

  // Compute matrix multiplication (mat1 * mat2) and return the result.
  // Return nullptr if dimensions mismatch for input matrices.
  static std::unique_ptr<RowMatrix<T>> MultiplyMatrices(std::unique_ptr<RowMatrix<T>> mat1,
                                                        std::unique_ptr<RowMatrix<T>> mat2) {
    // TODO(P0): Add code
    int row = mat1->GetRows();
    int c1 = mat1->GetColumns();
    int r2 = mat2->GetRows();
    int col = mat2->GetColumns();
    if (c1 != r2) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }
    auto mat_ptr = new RowMatrix<T>(row, col);

    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < col; ++j) {
        int sum = 0;
        for (int k = 0; k < c1; ++k) {
          sum += mat1->GetElem(i, k) * mat2->GetElem(k, j);
        }
        mat_ptr->SetElem(i, j, sum);
      }
    }

    return std::unique_ptr<RowMatrix<T>>(mat_ptr);
  }

  // Simplified GEMM (general matrix multiply) operation
  // Compute (matA * matB + matC). Return nullptr if dimensions mismatch for input matrices
  static std::unique_ptr<RowMatrix<T>> GemmMatrices(std::unique_ptr<RowMatrix<T>> matA,
                                                    std::unique_ptr<RowMatrix<T>> matB,
                                                    std::unique_ptr<RowMatrix<T>> matC) {
    // TODO(P0): Add code
    int ma_row = matA->GetRow();
    int ma_col = matB->GetColumns();
    int mb_row = matB->GetRow();
    int mb_col = matB->GetColumns();
    int mc_row = matC->GetRow();
    int mc_col = matC->GetColumns();
    if (ma_col != mb_row || ma_row != mc_row || mb_col != mc_col) {
      return std::unique_ptr<RowMatrix<T>>(nullptr);
    }
    auto mat_ptr = MultiplyMatrices(matA, matB);
    mat_ptr = AddMatrices(mat_ptr, matC);
    return std::unique_ptr<RowMatrix<T>>(mat_ptr);
  }
};
}  // namespace bustub
