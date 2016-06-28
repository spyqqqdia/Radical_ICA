package org.scalanlp.radical



import breeze.linalg.{
DenseMatrix, DenseVector, eigSym, svd,
max => bn_max, sum, Axis
}
import breeze.numerics.{
abs => bn_abs, sqrt => bn_sqrt
}



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



}


