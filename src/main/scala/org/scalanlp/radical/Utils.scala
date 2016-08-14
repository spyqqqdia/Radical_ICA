package org.scalanlp.radical



import breeze.linalg.{DenseMatrix, DenseVector, sum, max => bn_max}
import breeze.numerics.{abs => bn_abs}



object Utils {

   /**
     * @return Square matrix with diagonal d.
     */
    def diagMatrix(d:DenseVector[Double]):DenseMatrix[Double] = {

        val dim = d.length
        val A:DenseMatrix[Double] = DenseMatrix.zeros[Double](dim,dim)
        (0 until dim).map(i => A(i,i)=d(i) )
        A
    }

    /**
      * Reset A(i,j)=f(i,j).
      */
    def fillMatrix(A:DenseMatrix[Double],f:(Int,Int)=>Double):Unit =
        (0 until A.rows).map(i => (0 until A.cols).map(j => A(i,j)=f(i,j)))


    def time(execution: () => Unit) {
        val start = System.currentTimeMillis()
        execution()
        val time = System.currentTimeMillis() - start
        println("\nExecution time: " + time + "ms\n")
    }

    /**
      * Checks if A is a permutation matrix using the following criterion:
      * the maximal element of each row as well as the sum of all absolute values
      * of row entries differ from one by at most tol and the same is true of all
      * columns of A
      */
     def isPermutationMatrix(A:DenseMatrix[Double], tol:Double):Boolean = {

        var res = true
        (0 until A.rows).map(i => res &= bn_abs(bn_max(A(i,::))-1.0)<tol && bn_abs(sum(bn_abs(A(i,::)))-1.0)<tol)
        (0 until A.cols).map(j => res &= bn_abs(bn_max(A(::,j))-1.0)<tol && bn_abs(sum(bn_abs(A(::,j)))-1.0)<tol)
        res
    }
}


