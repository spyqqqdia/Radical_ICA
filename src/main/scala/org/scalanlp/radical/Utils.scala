package org.scalanlp.radical



import breeze.linalg.{DenseMatrix, DenseVector, sum, max => bn_max}
import breeze.numerics.abs



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
      * Checks if A is a permutation matrix.
      *
      * @param tol tolerated deviation of entries from 0 or 1.
      */
     def isPermutationMatrix(A:DenseMatrix[Double], tol:Double):Boolean = {

        var res = true
        val m = A.rows
        assert(m==A.cols,"A must be a square matrix")
        var a=0  // number of entries close to 1
        var b=0  // number of entries close to 0

        // check that each row has exaclty one 1.0, all other 0.0 entries
        for(i <- 0 until m){

            a=0; b=0
            for(j <- 0 until m) {
                if (abs(A(i,j) - 1.0) < tol) a += 1
                if (abs(A(i,j)) < tol) b += 1
            }
        }
        res &= (a==1) & (b==m-1)
        if(!res) return false;

        // now check the columns
        for(j <- 0 until m){

            a=0; b=0
            for(i <- 0 until m) {
                if (abs(A(i,j) - 1.0) < tol) a += 1
                if (abs(A(i,j)) < tol) b += 1
            }
        }
        res &= (a==1) & (b==m-1)
        res
    }
}


