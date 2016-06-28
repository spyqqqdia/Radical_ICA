package org.scalanlp.radical

import org.apache.xmlgraphics.ps.dsc.DefaultNestedDocumentHandler
import org.scalanlp.radical._

import breeze.linalg.{max => bn_max, _}
import breeze.numerics.{ abs => bn_abs}
import org.junit.Test
import junit.framework.TestCase
import org.junit.Assert._

/**
 * Created by oar on 4/10/16.
 */
class RotationTests extends TestCase {


    /** Allocate rotation with angles angle(0,1)=pi/4 and angle(1,2)=pi/4
      * and check that this rotates e_1 to e_3, e_2 to -e_1 and e_3 to -e_2
      */
    def testBasicRotation = {

        System.out.print("\nDoing RadicalTests.testBasicRotation: ")
        val pi = 3.1415926535897932
        val rot = new Rotation(3)

        rot.addRotation(0,1,pi/2)
        rot.addRotation(1,2,pi/2)

        val e1 = DenseVector[Double](1,0,0)
        val e2 = DenseVector[Double](0,1,0)
        val e3 = DenseVector[Double](0,0,1)

        val Q:DenseMatrix[Double] = rot.rotationMatrix

        val Qe1:DenseVector[Double] = Q*e1
        val Qe2:DenseVector[Double] = Q*e2
        val Qe3:DenseVector[Double] = Q*e3

        val diff = sum(bn_abs(Qe1-e3))+sum(bn_abs(Qe2+e1))+sum(bn_abs(Qe3+e2))
        assertTrue(diff < 1e-10)
        System.out.print("passed.\n")
    }

    /** Allocate a rotation and check if the inverse rotation matrix is computed correctly.
      * Rotation angles(i,j) are set to (i+j)/10.
      * Dimension of space in which we rotate is set to 20.
      */
    def testRotationInverse = {

        System.out.print("Doing RadicalTests.testRotationInverse: ")
        val dim=20
        val rot = new Rotation(dim)
        (0 until dim).map( i => ((i+1) until dim).map(j => {
            rot.addRotation(i,j,(i+j)/10.0)
        }))
        val M:DenseMatrix[Double] = rot.rotationMatrix
        val M_inv :DenseMatrix[Double] = rot.rotationMatrixInverse
        val Id:DenseMatrix[Double] = DenseMatrix.eye[Double](dim)

        val RQ:DenseMatrix[Double] = M*M_inv
        val QR:DenseMatrix[Double] = M_inv*M

        val leftDiff = bn_max(bn_abs(Id-RQ))
        val rightDiff = bn_max(bn_abs(Id-QR))

        assertTrue(leftDiff<1e-8 & rightDiff < 1e-8)
        System.out.print("passed.\n")
    }
    /** Allocate a rotation and check if the rotation matrix is orthogonal.
      * Rotation angles(i,j) are set to (i+j)/10.
      * Dimension of space in which we rotate is set to 20.
      */
    def testRotationIsOrthogonal = {

        System.out.print("Doing RadicalTests.testRotationIsOrthogonal: ")
        val dim = 20
        val rot = new Rotation(dim)
        (0 until dim).map( i => ((i+1) until dim).map(j => {
            rot.addRotation(i,j,(i+j)/10.0)
        }))
        val M:DenseMatrix[Double] = rot.rotationMatrix
        val Mt :DenseMatrix[Double] = M.t     // transpose
        val Id:DenseMatrix[Double] = DenseMatrix.eye[Double](dim)

        val RQ:DenseMatrix[Double] = M*Mt
        val QR:DenseMatrix[Double] = Mt*M

        val leftDiff = bn_max(bn_abs(Id-RQ))
        val rightDiff = bn_max(bn_abs(Id-QR))

        assertTrue(leftDiff<1e-8 & rightDiff < 1e-8)
        System.out.print("passed.\n")
    }
}
