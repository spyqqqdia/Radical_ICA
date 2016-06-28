package org.scalanlp.radical

import org.junit.Test
import junit.framework.TestCase
import org.junit.Assert._
import java.util.Random

import breeze.linalg.{DenseMatrix, DenseVector, eigSym, svd, max => bn_max}
import breeze.numerics.{abs => bn_abs}
import breeze.stats.distributions.Rand

/**
  * Created by oar on 6/3/16.
  */
class MatrixFunctionTests {

    val rng = new Random()
    val m=20
    val n=30
    // square matrix for decomposition tests
    val A = DenseMatrix.fill[Double](m,m){ rng.nextDouble() }
    // random data matrix, each col is a value of the random vector
    // correlation: multiplication by A
    val B = A*DenseMatrix.fill[Double](m,n){ rng.nextDouble() }


    // show upper left axb corner
    def showSelf(a:Int,b:Int) = {

        val msg = "\nMatrixFunctionTests: with "+m+"x"+m+" matrix A (upper left "+a+"x"+b+" corner):\n"
        println(msg)
        print(A(0 until a, 0 until b).toString+"\n\n")
    }
    /** check if A'A yields a symmetric matrix */
    def testSymmetry: Unit ={

        System.out.print("Testing symmetry of A*At: ")
        val S = A.t*A
        assertTrue(bn_max(bn_abs(S-S.t))<1e-10)    // ||S-St||_oo
        System.out.print("passed.\n")
    }
    /** eigenvector decomosition of symmetric matrix */
    def testSymEigVal: Unit = {

        System.out.print("Testing symmetric eigenvalue decomposition: ")
        val S = (A.t+A):*5.0
        // named parameters, must use "ev" (eigenvalues) and "rev" (right eigenvectors) here
        val eigSym.EigSym(ev,rev) = eigSym(S)
        val D = MatrixFunction.diag(ev)
        val Q =(rev*D)*rev.t

        assertTrue(bn_max(bn_abs(S-Q))<1e-10)   // ||S-Q||_oo
        System.out.print("passed.\n")
    }
    /** svd decomosition of A, warning: this has the form A=U*diag(s)*V not A=U*diag(s)*V' */
    def testSVD: Unit = {

        System.out.print("Testing SVD decmposition: ")
        val svd.SVD(u,s,v) = svd(A)
        val D = MatrixFunction.diag(s)
        val Q =(u*D)*v

        assertTrue(bn_max(bn_abs(A-Q))<1e-10)   // ||A-Q||_oo
        System.out.print("passed.\n")
    }
    /** check if the data whitening algo is correct */
    def testWhitening:Unit = {

        System.out.print("Testing data whitening: ")
        val D = MatrixFunction.whiten(B)
        val C = MatrixFunction.cov(D)
        // mxm identity matrix
        val I = DenseMatrix.eye[Double](m)

        assertTrue(bn_max(bn_abs(C-I))<1e-10)   // ||C-I||_oo
        System.out.print("passed.\n")
    }

    /** Allocates white and nearly white data and shows the covariance and whitening matrix.
      * Shows graphically that the whitening matrix is ony unique up to multiplication with a
      * orthogonal matrix. Even if we start with very close to white data the whitening matrix need
      * not be close to the identity.
      *
      * If the covariance matrix is exactly the identity, then so is the whitening matrix.
      * But if the covariance matrix is differs slightly from the identity, then the whitening matrix
      * can deviate significantly from the identity (but will be close to orthogonal).
      */
    def whiteningTest(): Unit ={

        val dim = 3
        System.out.println("\nFirst nearly white data:")
        var X: DenseMatrix[Double] = Math.sqrt(3.0)*(2.0*DenseMatrix.rand[Double](dim,1000,Rand.uniform)-1.0)
        var C = MatrixFunction.cov(X)
        System.out.println("Covariance matrix:\n"+C+"\n")
        var Q = MatrixFunction.whiteningMatrix(X)
        System.out.println("\nWhitening matrix:\n"+Q+"\n")
        System.out.println("\nCovariance matrix after whitening:\n"+MatrixFunction.cov(Q*X)+"\n")

        System.out.println("\n\nFully white data:")
        X = MatrixFunction.whiten(2.0*DenseMatrix.rand[Double](dim,1000,Rand.uniform)-1.0)
        C = MatrixFunction.cov(X)
        System.out.println("Covariance matrix:\n"+C+"\n")
        Q = MatrixFunction.whiteningMatrix(X)
        System.out.println("\nWhitening matrix:\n"+Q+"\n")
        System.out.println("\nCovariance matrix after whitening:\n"+MatrixFunction.cov(Q*X)+"\n")
    }

}
