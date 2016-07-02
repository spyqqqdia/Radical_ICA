package org.scalanlp.radical

import scala.collection.mutable.ListBuffer
import breeze.linalg.{Axis, DenseMatrix, DenseVector, eigSym, sum, svd, max => bn_max}
import breeze.numerics.{abs => bn_abs, sqrt => bn_sqrt}
import Utils._
import org.junit.Assert._



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
        // prepend to get them in reverse order for computation of the inverse
        ((i,j,phi)) +=: jacobiRotations
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

    /** Inverse of the rotation matrix, alternative computation. Because of orthogonality
      * the inverse is simply the transpose so this method is redundant. It's useful for consistency checks
      * though.
      */
    def rotationMatrixInverse:DenseMatrix[Double]= {

        val rotInv = new Rotation(dim)
        // this works only since jacobiRotations contains them in reverse order
        jacobiRotations.map(t => rotInv.addRotation(t._1,t._2,-t._3))
        rotInv.rotationMatrix
    }

   /**
    * Checks if the rotation matrix is orthogonal and if the inverse matrix is computed correctly.
    */
    def selfTest:Boolean = {

        val I = DenseMatrix.eye[Double](dim)
        val R =rotationMatrix
        val Q :DenseMatrix[Double] = R.t     // transpose

        var leftDiff = bn_max(bn_abs(I-R*Q))
        var rightDiff = bn_max(bn_abs(I-Q*R))
        assert(
            leftDiff<1e-8 & rightDiff < 1e-8,
            "Rotation matrix is not orthogonal."
        )
        val R_inv = rotationMatrixInverse
        leftDiff = bn_max(bn_abs(I-R*R_inv))
        rightDiff = bn_max(bn_abs(I-R_inv*R))
        assert(
            leftDiff<1e-8 & rightDiff < 1e-8,
            "Inverse of rotation matrix not correct."
        )
        print("passed.\n")
        true
    }
    /** Note that J(i,j,alpha)*J(i,j,beta)=J(i,j,alpha+beta). This method reduces the list of Jacobi rotations
      * in this rotation applying this equation until each coordinate pair (i,j) occurs only one.
      * No guarantees are given about the order in which the resulting Jacobi rotations are applied.
      *
      * WARNING: in general the resulting rotation is only guaranteed to be equivalent to this rotation
      * (same rotation matrix) if all the Jacobi rotations in this rotation _commute_ pairwise.
      */
    def reduced:Rotation = {

        val reducedJacobiRotations:Iterable[(Int,Int,Double)] = jacobiRotations.groupBy(t=>(t._1,t._2)).map({

            case ((i,j),rotationList) => (i,j,rotationList.map(J=>J._3).sum)
        })
        val rot = new Rotation(dim)
        val it = reducedJacobiRotations.iterator
        while(it.hasNext){ val J = it.next; rot.addRotation(J._1,J._2,J._3); }
        rot
    }
}





object Rotation {

    /** Factory method to allocate rotation object.*/
    def apply(dim:Integer):Rotation = new Rotation(dim)
}
