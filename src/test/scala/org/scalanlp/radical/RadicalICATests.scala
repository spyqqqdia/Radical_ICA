package org.scalanlp.radical



import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions.{ContinuousDistr,Uniform}
import scala.collection.mutable.ListBuffer

/**
  * Created by oar on 6/4/16.
  */
class RadicalICATests {

    val pi = 3.1415926535897932
    val rng = scala.util.Random

    /** Pad 3x2 matrix with all entries = 1.0 with nPad=3 and sigmaPad=0.1.
      * Then print result. Must be 3x6 matrix with entries close to 1.0.
      */
    def showPadding = {

        val data = DenseMatrix.fill[Double](5,2){1.0};
        val paddedData = RadicalICA.pad(data,3,0.1);
        System.out.println(MathTools.applyFcn(paddedData,(u:Double)=>MathTools.round(u,2)));
    }
    /**
      * The test allocates a data sample of the type [[DataGenerator.rotatedWhitenedSample]] with
      * rotation rot rotating the sample by -0.7 radians in non interacting dimensions (0,1),(2,3),...(d-2,d-1).
      * Here d=dim=coordinateDistributions.size is the data dimension. Note that this rotation matrix is rather
      * sparse.
      *
      * With this the demixing rotation is uniquely determined  and must rotate by 0.7 radians in the same
      * dimensions (0,1), (2,3),...,(d-2,d-1).
      * Radical ICA is now carried out and we check if the reverse rotations are all discovered. Note that the
      * rotation matrix is very sparse.
      *
      * In verbose mode we print all the Jacobi rotations of the demixing
      * rotation with angle |phi| > precision so that we can check how close we are to this solution.
      *
      * Because of sampling error (distribution of the sample deviates from the true distribution generating the sample)
      * it is not easy to give a criterion for failing the test. Here we are reporting failure if the demixing rotation
      * contains Jacobi rotations with angle phi off target by more than precision radians. You can set the precision
      * as a parameter.
      *
      * The failure threshold depends on the data dimension, the coordinate distributions and the sample size,
      * but should be set significantly smaller than the size 0.7 of the actual Jacobi rotations.
      *
      * With this the test passes if and only if the demixing rotation contains all the Jacobi rotations
      * J(2i,2i+1,phi) with |phi-0.7| < precision and all other Jacobi rotations which it contains have angle
      * |phi| < precision.
      *
      * @param coordinateDistributions list of dim distributions assigned to components of random vector,
      *                                size dim must be even.
      * @param sampleSize size of padded sample (unpadded size will be sampleSize/40).
      * @param entropyEstimator  estimator for the marginal entropies. Provided are
      *                          MathTools.entropyEstimatorVasicek and MathTools.entropyEstimatorEmpirical.
      * @param doParallelSearch  use the multithreaded version of the rotation search,
      *                          the sequential version exists only for testing and debugging.
      * @param precision  acceptable deviation (in radians) of angles in Jacobi rotations from theoretical optimum.
      * @param verbose  if set to true will show all Jacobi rotations J(i,j,phi) of the demixing matrix where
      *                 |phi| > precision, otherwise will only report pass/fail as defined above.
      *
      *
      * */
    def testRadicalSparseWithEstimator(
            coordinateDistributions:List[ContinuousDistr[Double]], sampleSize:Int,
            entropyEstimator:(DenseVector[Double])=>Double,
            doParallelSearch:Boolean, precision:Double, verbose: Boolean
    ): Boolean = {

        val dim = coordinateDistributions.size
        assert(dim%2==0,"Data dimension (coordinateDistributions.size) must be even but is = "+dim)

        val rot = new Rotation(dim)
        // rotate in noninteracting coordinates, so the inverse rotation is unique
        (0 until dim/2).map(i => rot.addRotation(2*i,2*i+1,-0.7))

        val sigma = 0.15
        val data = DataGenerator.rotatedWhitenedSample(coordinateDistributions,sampleSize,sigma,rot)

        val nAngles = 100
        // important: no whitening or else the whitening matrix messes up our rotation.
        val doWhitenData = false

        Timer start
        val ica = RadicalICA(data,nAngles,entropyEstimator,doWhitenData,doParallelSearch,verbose)
        val W: Rotation = ica.rot
        Timer.stop
        System.out.println(Timer.report)

        // since the jacobi rotations in this example commute pairwise, we can reduce the rotation
        // to eliminate multiple rotations in the same coordinate plane (i,j)
        val V:Rotation = W.reduced
        // go through the list of Jacobi rotations of W to asses pass/fail
        val jacobiRotations:ListBuffer[(Int,Int,Double)]  = V.jacobiRotations

        // list of test results true/false
        val results = jacobiRotations.map(
            t => if(t._1%2==0 & t._2==1+t._1) Math.abs(t._3-0.7)<=precision else Math.abs(t._3)<=precision
        )
        // did all of them pass
        val pass = results.foldLeft(true)((a:Boolean,b:Boolean)=>a&b)

        println((if(pass) "Test passed" else "Test failed")+" with precision "+precision+".")
        if(verbose){
            System.out.println("\nOptimal rotation found:")
            System.out.print(W.listRotationAngles(precision))
        }
        pass
    }
    /** Does the test [[testRadicalSparseWithEstimator]] for
      * entropyEstimator = [[MathTools.entropyEstimatorEmpiricalFast]], [[MathTools.entropyEstimatorEmpiricalAdapted]]
      * and [[MathTools.entropyEstimatorVasicek]] where the distribution is Uniform(-a,a) with a=sqrt(3) (thus standard
      * deviation = 1.0) in each coordinate.
      *
      * This is mainly a test of the speed of the algorithm with the various estimators. For a test with more
      * interesting distributions see [[testRadicalSparse]].
      *
      * @param sampleSize size of padded sample (unpadded size will be sampleSize/40).
      * @param doVasicek if set to true the test will also be run with the Vasicek entropy estimator.
      * In high dimensions (10+) at the usual sample sizes (10000+) this is very slow. The test always runs with
      * the fast and adapted empirical entropy estimator.
      * @param doParallelSearch use the multithreaded version.
      * @param precision acceptable deviation (in radians) of angles in Jacobi rotations from theoretical optimum.
      * @param verbose if set to true prints the Jacobi angles of the demixing matrix.
      */
    def testCube(dim:Int, sampleSize:Int, doParallelSearch:Boolean, doVasicek:Boolean, precision:Double, verbose: Boolean) = {

        var msg = "\nData X: uniform in d-dimensional cube [-a,a]^d with d="+dim+" and a chosen such that cov(X)=Id.\n"
        msg += "Rotating data with Jacobi rotations J(0,1,-0.7), J(2,3,-0.7),..., J(d-2,d-1,-0.7).\n"
        msg += "Note: coordinates are noninteracting with respect to these rotations.\n"
        msg += "Thus the unique product of Jacobi rotations J(i,j,phi) with angles phi in [-pi/4,pi/4]\n"
        msg += "undoing the data rotation consists of the inverse of these Jacobi rotations\n"
        msg += "(sign of the angle reversed).\n"
        msg += "Our ICA algorithm must find just this rotation.\n"
        if(verbose) {
            msg += "We will print the rotations J(i,j,phi) found, showing only those for which "
            msg += "|phi| > "+precision+" (radians)."
        }
        System.out.println(msg+"\n\n")

        val a = Math.sqrt(3)
        val dists = ListBuffer[ContinuousDistr[Double]]()
        for(i <- 0 until dim) dists += Uniform(-a,a)
        val coordinateDistributions = dists.toList

        System.out.println("\n####-----Cube rotation test with empirical entropy estimator-----####\n")
        var entropyEstimator = MathTools.entropyEstimatorEmpiricalFast
        testRadicalSparseWithEstimator(
            coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose
        )


        System.out.println("\n####-----Cube rotation test with adapted empirical entropy estimator-----####\n")
        entropyEstimator = MathTools.entropyEstimatorEmpiricalAdapted(sampleSize/50)
        testRadicalSparseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose)

        if(doVasicek){

            System.out.println("\n####-----Cube rotation test with Vasicek entropy estimator-----####\n")
            val m = 200  // spacing in the Vasicek estimator appropriate for the samples in testCubeWithEstimator
            entropyEstimator = MathTools.entropyEstimatorVasicek(m)
            testRadicalSparseWithEstimator(
                coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose
            )

        }
    }
    /** Does the test [[testRadicalSparseWithEstimator]] for
      * entropyEstimator = [[MathTools.entropyEstimatorEmpiricalFast]], [[MathTools.entropyEstimatorEmpiricalAdapted]]
      * and [[MathTools.entropyEstimatorVasicek]].
      *
      * @param sampleSize size of padded sample (unpadded size will be sampleSize/40).
      * @param doVasicek if set to true the test will also be run with the Vasicek entropy estimator.
      * In high dimensions (10+) at the usual sample sizes (10000+) this is very slow. The test always runs with
      * the fast and adapted empirical entropy estimator.
      * @param doParallelSearch use the multithreaded version.
      * @param precision acceptable deviation (in radians) of angles in Jacobi rotations from theoretical optimum.
      * @param verbose if set to true prints the Jacobi angles of the demixing matrix.
      */
    def testRadicalSparse(
           coordinateDistributions:List[ContinuousDistr[Double]], sampleSize:Int, doParallelSearch:Boolean,
           doVasicek:Boolean, precision:Double, verbose: Boolean
     ) = {

        var msg = "\nData X with independent coordinates and assigned distribution in each coordinate\n"
        msg += "are whitended and rotated with Jacobi rotations J(0,1,-0.7), J(2,3,-0.7),..., J(d-2,d-1,-0.7).\n"
        msg += "Note: coordinates are noninteracting with respect to these rotations.\n"
        msg += "Thus the unique product of Jacobi rotations J(i,j,phi) with angles phi in [-pi/4,pi/4]\n"
        msg += "undoing the data rotation consists of the inverse of these Jacobi rotations\n"
        msg += "(sign of the angle reversed).\n"
        msg += "Our ICA algorithm must find just this rotation.\n"
        if(verbose) {
            msg += "We will print the rotations J(i,j,phi) found, showing only those for which "
            msg += "|phi| > "+precision+" (radians)."
        }
        System.out.println(msg+"\n")

        System.out.println("\n\n####-----Cube sparse rotation test with empirical entropy estimator-----####")
        var entropyEstimator = MathTools.entropyEstimatorEmpiricalFast
        testRadicalSparseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose)


        System.out.println("\n\n####-----Cube sparse rotation test with adapted empirical entropy estimator-----####")
        entropyEstimator = MathTools.entropyEstimatorEmpiricalAdapted(sampleSize/50)
        testRadicalSparseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose)

        if(doVasicek){

            System.out.println("\n\n####-----Cube sparse rotation test with Vasicek entropy estimator-----####")
            val m = 200  // spacing in the Vasicek estimator appropriate for the samples in testCubeWithEstimator
            entropyEstimator = MathTools.entropyEstimatorVasicek(m)
            testRadicalSparseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,precision,verbose)

        }
    }

    /**
      * The test allocates a data sample of the type [[DataGenerator.rotatedWhitenedSample]] with a rotation
      * affecting all coordinate pairs, each with a random angle u*pi, where u in (0,1) is random. Thus the rotation
      * matrix U is dense.
      *
      * Because the data are white the demixing matrix coincides with the rotation matrix V computed by the Radical
      * algorithm. It is not uniquely determined but the product VU with the rotation matrix U applied to
      * the independent data sample with  must be a permutation matrix. This condition is checked.
      *
      * Note that VU is the transformation which first rotates the independent sources with U, then reverses this
      * rotation with the rotation V found by the Radical algorithm. If the algorithm works perfectly this
      * transformation can only permute coordinates. The same is not true of the transformation UV, i.e. the
      * matrix UV need not be a permutation matrix.
      *
      * The algorithm is not expected to work in large dimensions if the independent source are rotated with
      * a dense rotation matrix. Thus keep the dimension low. In verbose mode we print the product VU.
      *
      * @param coordinateDistributions list of dim distributions assigned to components of random vector in order,
      *      see [[DataGenerator.rotatedWhitenedSample]]. The size of this list is the data dimension.
      * @param sampleSize size of padded sample (unpadded size will be sampleSize/40).
      * @param entropyEstimator  estimator for the marginal entropies. Provided are
      *      MathTools.entropyEstimatorVasicek and MathTools.entropyEstimatorEmpirical.
      * @param doParallelSearch  use the multithreaded version of the rotation search,
      *      the sequential version exists only for testing and debugging.
      * @param tolerance  acceptable deviation from permutation matrix, see [[Utils.isPermutationMatrix]].
      * @param verbose  if set to true will print the product UV which should be a permutation matrix.
      *
      *
      * */
    def testRadicalDenseWithEstimator(
            coordinateDistributions:List[ContinuousDistr[Double]], sampleSize:Int,
            entropyEstimator:(DenseVector[Double])=>Double,
            doParallelSearch:Boolean, tolerance:Double, verbose: Boolean
        ): Boolean = {

        val dim = coordinateDistributions.size
        assert(dim%2==0,"Data dimension (coordinateDistributions.size) must be even but is = "+dim)

        val rot = new Rotation(dim)
        // rotate in noninteracting coordinates, so the inverse rotation is unique
        for(i <- 0 until dim; j <- (i+1) until dim) rot.addRotation(i,j,pi*rng.nextFloat())
        val U:DenseMatrix[Double] = rot.rotationMatrix

        val sigma = 0.15
        val data = DataGenerator.rotatedWhitenedSample(coordinateDistributions,sampleSize,sigma,rot)

        val nAngles = 100
        // important: no whitening or else the whitening matrix messes up our rotation.
        val doWhitenData = false

        Timer start
        val ica = RadicalICA(data,nAngles,entropyEstimator,doWhitenData,doParallelSearch,verbose)
        val V: DenseMatrix[Double] = ica.rot.rotationMatrix
        Timer.stop
        System.out.println(Timer.report)

        val VU:DenseMatrix[Double] = V*U;

        // did all of them pass
        val pass = Utils.isPermutationMatrix(VU,tolerance)

        println((if(pass) "Test passed" else "Test failed")+" with precision "+tolerance+".")
        if(verbose){
            System.out.println("\nMatrix UV, should be a permuation matrix:")
            System.out.print(MathTools.round(VU,3))
            assert(MathTools.isOrthogonal(VU,1e-13))
        }
        pass
    }
    /** Does the test [[testRadicalDenseWithEstimator]] for
      * entropyEstimator = [[MathTools.entropyEstimatorEmpiricalFast]], [[MathTools.entropyEstimatorEmpiricalAdapted]]
      * and [[MathTools.entropyEstimatorVasicek]].
      *
      * @param coordinateDistributions list of dim distributions assigned to components of random vector in order,
      *      see [[DataGenerator.rotatedWhitenedSample]]. The size of this list is the data dimension.
      * @param sampleSize size of padded sample (unpadded size will be sampleSize/40).
      * @param doVasicek if set to true the test will also be run with the Vasicek entropy estimator.
      * In high dimensions (10+) at the usual sample sizes (10000+) this is very slow. The test always runs with
      * the fast and adapted empirical entropy estimator.
      * @param doParallelSearch use the multithreaded version.
      * @param tolerance acceptable deviation from permutation matrix, see see [[Utils.isPermutationMatrix]].
      * @param verbose if set to true prints the Jacobi angles of the demixing matrix.
      */
    def testRadicalDense(
            coordinateDistributions:List[ContinuousDistr[Double]], sampleSize:Int, doParallelSearch:Boolean,
            doVasicek:Boolean, tolerance:Double, verbose: Boolean
        ) = {

        var msg = "\nData X with independent coordinates and assigned distribution in each coordinate\n"
        msg += "are whitended and rotated with dense rotation matrix U.\n"
        msg += "If V is the demixing rotation computed by Radical, then UV must be a permutation matrix.\n"
        msg += "Will check this now."
        System.out.println(msg+"\n")

        System.out.println("\n\n####-----Cube dense rotation test with empirical entropy estimator-----####")
        var entropyEstimator = MathTools.entropyEstimatorEmpiricalFast
        testRadicalDenseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,tolerance,verbose)


        System.out.println("\n\n####-----Cube dense rotation test with adapted empirical entropy estimator-----####")
        entropyEstimator = MathTools.entropyEstimatorEmpiricalAdapted(sampleSize/50)
        testRadicalDenseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,tolerance,verbose)

        if(doVasicek){

            System.out.println("\n\n####-----Cube dense rotation test with Vasicek entropy estimator-----####")
            val m = 200  // spacing in the Vasicek estimator appropriate for the samples in testCubeWithEstimator
            entropyEstimator = MathTools.entropyEstimatorVasicek(m)
            testRadicalDenseWithEstimator(coordinateDistributions,sampleSize,entropyEstimator,doParallelSearch,tolerance,verbose)
        }
    }


}
