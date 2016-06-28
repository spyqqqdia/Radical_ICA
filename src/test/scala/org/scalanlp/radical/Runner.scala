package org.scalanlp.radical


/**
 * Created by oar on 4/10/16.
 * Contains main method to run the various tests.
 */
object Runner {

    val rotationTests = new RotationTests
    val matrixFunctionTests = new MatrixFunctionTests
    val radICATests = new RadicalICATests

    val runTests = false


    def main (args: Array[String]) {

        print("\nRunner: starting...")

        // unit tests
        if (runTests) {
            rotationTests.testBasicRotation
            rotationTests.testRotationConsistency

            matrixFunctionTests.testSymmetry
            matrixFunctionTests.testSymEigVal
            matrixFunctionTests.testSVD
            matrixFunctionTests.testWhitening
        }

        // ad hoc tests
        // radICATests.showPadding
        // matrixFunctionTests.whiteningTest()
        // MathTests.entropyTest()


        val dim = 10
        val precision = 0.08
        val doVasicek = false
        val doParallelSearch =true
        val verbose = true
        radICATests.testCube(dim,doParallelSearch,precision,doVasicek,verbose)


    }
}
