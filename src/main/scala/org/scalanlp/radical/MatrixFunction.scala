package org.scalanlp.radical

import breeze.linalg.{Axis, DenseMatrix, DenseVector, eigSym, sum }
import breeze.numerics.{ sqrt => bn_sqrt }


/**
  * Created by oar on 6/3/16.
  */
/** Some static helper functions for matrices
  *
  */
object MatrixFunction {

    /**
      * @param d diagonal of matrix
      * @return diagonal matrix with diagonal d
      */
    def diag(d:DenseVector[Double]): DenseMatrix[Double] = {

        val D = DenseMatrix.zeros[Double](d.length,d.length)
        (0 until d.length).map(i => D(i,i)=d(i))
        D
    }
    /**
      * @param A data matrix, each column is a value of the random vector,
      *          each row is a sample of a (coordinate) random variable,
      * @return  empiricial covariance matrix of data
      */
    def cov(A:DenseMatrix[Double]):DenseMatrix[Double] = {

        // set col means to zero
        val n = A.cols
        val D:DenseMatrix[Double] = A.copy
        val mu:DenseVector[Double] = sum(D,Axis._1):*(1.0/n)    // sum along rows --> col vector
        (0 until n).map(i => D(::,i):-=mu)
        val C = (D*D.t):*(1.0/n)
        // make exactly symmetric
        (C+C.t):*(0.5)
    }

    /**
      * @param A data matrix, each column is a value of the random vector,
      *          each row is a sample of a (coordinate) random variable,
      * @return  matrix Q such that cov(QA)=I computed from eigenvalue decomposition.
      *          Note that Q is determined only up to multiplication by an orthogonal matrix.
      */
    def whiteningMatrix(A:DenseMatrix[Double]):DenseMatrix[Double] = {

        val eigSym.EigSym(l,u) = eigSym(cov(A))          // l: eigenvalues, u: right eigenvectors

        // now  CovM = U*diag(l)*U'  and so we want I = cov(QA) = Qcov(A)Q' = (QU)*diag(l)*(QU)'
        // i.e. QU=1/sqrt(diag(l)))   i.e   Q = diag(1/sqrt(l))U'
        val w = DenseVector.ones[Double](l.length):/bn_sqrt(l)
        diag(w)*u.t

    }
    /**
      * @param A data matrix, each column is a value of the random vector,
      *          each row is a sample of a (coordinate) random variable,
      * @return  matrix A made white: QA, where cov(QA)=I
      */
    def whiten(A:DenseMatrix[Double]):DenseMatrix[Double] = whiteningMatrix(A)*A
}