package org.scalanlp.radical

import scala.collection.mutable.ListBuffer
import breeze.linalg.{
DenseMatrix, DenseVector, eigSym, svd,
max => bn_max, sum, Axis
}
import breeze.numerics.{
abs => bn_abs, sqrt => bn_sqrt
}
import Utils._



/** Class to generate Rotations (orthogonal matrix with determinant one) as a product of Jacoby rotations
  *  J(i,j,phi).
  *
  *  J(i,j,phi) rotates counterclockwise by the angle phi (radians) in the two dimensional subspace
  *  spanned by the standard unit vectors e_i, e_j (i.e.: in coordinates i,j).
  *
  *  Such a rotation affects only coordinates i and j of any vector.
  *  Note that two such rotations commute in general only if they act on disjoint coordinates.
  *
  */
class Rotation(val dim:Int) {

    assert(dim>=1)

    /** List of rotations J(i,j,alpha) in the _reverse_ _order_ of application.
      * In this order we can easily compute the inverse of the rotation matrix.
      */
    val jacobiRotations:ListBuffer[(Int,Int,Double)] = new ListBuffer[(Int,Int,Double)]
    val rotationMatrix:DenseMatrix[Double]= DenseMatrix.eye[Double](dim)
    /** Add rotation by angle phi in coordinates (i,j), must satisfy 0<=i,j<dim and i != j.
     */
    def addRotation(i:Int,j:Int,phi:Double):Unit = {

        assert(0<=i && 0<=j && i<dim && j<dim && i!=j)
        ((i,j,phi)) +=: jacobiRotations     // prepend to get them in reverse order
        val a = Math.cos(phi)
        val b = Math.sin(phi)

        val row_i = DenseVector.zeros[Double](dim)
        val row_j = DenseVector.zeros[Double](dim)
        (0 until dim).map(k => {
            row_i(k) = rotationMatrix(i,k)
            row_j(k) = rotationMatrix(j,k)
        })
        rotationMatrix(i,::) := (a*row_i - b*row_j).t
        rotationMatrix(j,::) := (b*row_i + a*row_j).t
    }
    /** The orthogonal matrix implementing the rotation by multiplication on the left.
      */
    def apply():DenseMatrix[Double] = rotationMatrix
    /** Reset rotation to identity matrix.
      */
    def clear():Unit = {

      jacobiRotations.dropWhile(el=>true)    // clear
      Utils.fillMatrix(rotationMatrix, (i: Int, j: Int) => if (i == j) 1.0 else 0.0)
    }

    /**
      * @param cutoff only Jacobi rotations J(i,j,phi) with |phi|>cutoff are shown
      * @return String showing all the Jacobi rotations J(i,j,phi) as "angle(i,j)=phi".
     */
    def listRotationAngles(cutoff:Double=0.0):String = {

        val rotations = jacobiRotations.filter(t=>Math.abs(t._3)>cutoff)
        val vs:Vector[String] = rotations.map(t => "angle("+t._1+","+t._2+")="+MathTools.round(t._3,4)+"\n").toVector
        vs.foldLeft("\n")((s1:String,s2:String)=>s1+s2)
    }
    def rotationMatrixInverse:DenseMatrix[Double]= {

        val rotInv = new Rotation(dim)
        jacobiRotations.map(t => rotInv.addRotation(t._1,t._2,-t._3))
        rotInv.rotationMatrix
    }

   /**
    * Checks if the rotation matrix is orthogoanl
    */
    def selfTest:Boolean = {

        val I = DenseMatrix.eye[Double](dim)
        apply() == I
    }
}





object Rotation {

    /** Factory method to allocate rotation object.*/
    def apply(dim:Integer):Rotation = new Rotation(dim)
}
