package org.scalanlp.radical

import breeze.linalg.{Axis, DenseMatrix, DenseVector, eigSym, sum, svd, max => bn_max}
import scala.collection.mutable.ListBuffer
import scala.concurrent.{Await, Future}
import scala.concurrent.duration._
import scala.concurrent.ExecutionContext.Implicits.global



/**
  * Created by oar on 6/3/16.
  * Class to perform the ICA analysis on a data sample X as described in docs/RADICAL_ICA.pdf
  *
  * The model for the data X is as follows: X=AS, where S is a random vector with independent components
  * and A an invertible linear transformation (i.e. matrix). The components S_i of S are called the "sources"
  * and A is called the mixing matrix. You are interested in the demixing matrix B=inv(A) with which you recover
  * the sources S from the data X as S=BX.
  *
  * The algorithm proceeds as follows: the data are whitened and then padded:
  * each item in the whitened sample is multiplied nPad times and then a small random perturbation
  * drawn independently from N(0,sigma²) in each coordinate is added to it.
  *
  * This means that the sample distribution is convolved with the multinormal distribution N(0,sigma²I),
  * where I is the identity matrix of the same dimension as the data.
  * This smooths over sample idiosyncracies which can lead to spurious minima in the contrast function,
  * see docs/RADICAL_ICA.pdf, p3, section 3).
  *
  * Let pX denote the padded data and Y=QpX the whitened data sample, where Q is an invertible matrix such that
  * cov(Y)=I, the identity matrix (the augmented data matrix XX is of full rank with probability one).
  *
  * The algorithm now looks for a _rotation_ W such that the vector Z=WY has independent components.
  * To find such W it minimizes a contrast function F=F(W) which is the Kullback-Leibler distance of
  * the empirical distribution of Z from the product of the marginal distributions (distance from independence).
  * It can be shown that
  *
  *    F(Y) = const + sum_i entropy(Y_i),
  *
  * where Y_i denotes the i-th coordinate of the random vector Y.
  * Thus the algorithm tries to find the rotation W such that Y=WX minimizes the sum of the entropies
  * H_i=H(Y_i)=entropy(Y_i).
  * Note that the Y_i-sample is the the i-th row of the data matrix Y.
  *
  * The algorithm needs samples of large size (10000+). The usual implementation automatically pads the data and then
  * whitens the resulting sample. This procedure interferes with some of our tests.
  *
  * The tests rotate a set of independent signals by a given rotation U and check if the rotation V computed by the
  * algorithm reverses the rotation U. In general one can only guarantee VU must be permutation matrix.
  * If U is a product of Jacobi rotations J in non interacting (disjoint) coordinate pairs (i,j) (and thus necessarily
  * fairly sparse), then V has to reverse the rotations J one by one and some tests check this behaviour.
  *
  * For some of the tests to work properly the rotation U has to be applied after padding and whitening.
  * For this reason we include a switch to turn off whitening. It is then assumed that the data are white.
  * Likewise it is left to you to pad the data sample to reasonable size using the function
  * DataGenerator.smoothMultiDimensionalData. In this fashion we can run the tests on large white test data samples
  * which are not padded.
  *
  * @param data sample of multidimensional distribution, each column is one random value of the
  *              distribution, each row a sample of a coordinate variable. Sample size needs to be 10000+,
  *              use [[DataGenerator.smoothMultiDimensionalData]] to increase the sample size if needed.
  * @param nAngles number of equidistant grid points in [-pi/4,pi/4] scanned for optimal angles alpha of the
  *                Jacobi rotations (100 should be enough).
  * @param entropyEstimator estimator for the marginal entropies. Provided are
  *      [[MathTools.entropyEstimatorVasicek]], [[MathTools.entropyEstimatorEmpiricalAdapted]] and
  *      [[MathTools.entropyEstimatorEmpiricalFast]]
  * listed from slowest and most precise to fastest and least precise. The last estimator can fail badly on
  * strongly clustered data (spiky densities).
  *
  * @param doWhitenData padded data will be whitened before proceeding. Should always be done.
  *      But for testing purposes we want to be able to switch this off
  * @param doParallelSearch use the multithreaded version of the rotation search,
  *      the sequential version exists only for testing and debugging.
  */
class RadicalICA(
            val data:DenseMatrix[Double],
            val nAngles:Integer,
            val entropyEstimator:(DenseVector[Double])=>Double,
            val doWhitenData:Boolean = true,
            val doParallelSearch:Boolean = true,
            val verbose:Boolean = false
){
    val pi = 3.1415926535897932
    val dim = data.rows
    val Q = if(doWhitenData) MatrixFunction.whiteningMatrix(data) else DenseMatrix.eye[Double](dim)
    val X = if(doWhitenData) Q*data else data

    /** Rotation to make the whitened data independent.*/
    val rot:Rotation = if(doParallelSearch)
                            findOptimalRotationPar(verbose)
                       else
                            findOptimalRotation(verbose)
    /** Matrix B such that S=B*X is the sample of the independent components: since
      * Z=WY=WQX has independent components, we can take B = WQ.
      */
    val demixingMatrix:DenseMatrix[Double] = rot.rotationMatrix*Q

    /**
      * Performs complete grid search over [[nAngles]] equidistant angles in [-pi/4,pi/4] to find the optimal angle phi
      * for the Jacobi rotation J(i,j,phi) for each pair of coordinates (i,j) and  produces the demixing matrix
      * as product of the optimal Jacobi rotations. Sequential version.
      */
    private def findOptimalRotation(vrbose:Boolean):Rotation = {

        if(vrbose) System.out.println("\nRadicalICA: computing optimal rotation.")
        var rot = new Rotation(dim)

        // entropies = entropy(row_i(X))
        val entropies = new DenseVector[Double](dim)
        (0 until dim).map(i => entropies(i) = entropyEstimator(X(i,::).t))


        // at most 10 sweeps with Jacoby rotations J(i,j,alpha) over all coordinate pairs (i,j) with i<j
        // target: maximum entropy = minimum contrast function
        var sweep = 0
        var contrastFcn = sum(entropies)
        if(verbose) System.out.print("\nContrast function at start: "+MathTools.round(contrastFcn,2)+"\n")
        var contrastFcn_new = contrastFcn
        var delta_CF = Double.MinValue    // to enter the loop
        while(sweep<10 & delta_CF < -Math.abs(0.01*contrastFcn)){

            var i=0
            while(i<dim){
                var j=i+1
                val row_i:DenseVector[Double] = X(i,::).t   // row_i(X)
                while(j<dim){

                    val row_j:DenseVector[Double] = X(j,::).t   // row_j(X)
                    var dPhi=pi/(2*nAngles)
                    val (phi:Double,dCF:Double) = optimalRotation(i,j,-pi/4,dPhi,X,entropies)

                    if(dCF<0) {
                        // optimal angle and entropy decrease for coordinates (i,j) are set,
                        // apply rotation and update state to current optimum
                        val a = Math.cos(phi)
                        val b = Math.sin(phi)
                        val u_i = a*row_i-b*row_j
                        val u_j = b*row_i+a*row_j
                        X(i,::) := u_i.t
                        X(j,::) := u_j.t
                        entropies(i) = entropyEstimator(u_i)
                        entropies(j) = entropyEstimator(u_j)
                        rot.addRotation(i,j,phi)
                    }
                    j+=1
                } // while j
                i+=1
            } // while i
            contrastFcn_new = sum(entropies)
            delta_CF = contrastFcn_new-contrastFcn
            if(vrbose){
                val distDecreasePcnt = MathTools.round(100*delta_CF/contrastFcn,2)
                var msg:String = "Sweep "+sweep+": contrast function = "+MathTools.round(contrastFcn_new,2)
                msg += ", decrease (%) = "+distDecreasePcnt
                System.out.println(msg)
            }
            contrastFcn = contrastFcn_new
            sweep+=1
        } // while sweep
        rot
    } // findOptimalRotation


    /**
      * Performs complete grid search over [[nAngles]] equidistant angles in [-pi/4,pi/4] to find the optimal
      * angle phi for the Jacobi rotation J(i,j,phi) for each pair of coordinates (i,j) and  produces the demixing
      * matrix as product of the optimal Jacobi rotations, parallelized version.
      */
    private def findOptimalRotationPar(vrbose:Boolean):Rotation = {

        if(vrbose) System.out.println("\nRadicalICA: computing optimal rotation, parallelized.")
        val dPhi = pi/(2*nAngles)           // step width phi_step in [-pi/4,pi/4]
        var rot = new Rotation(dim)

        // entropies = entropy(row_i(X))
        val entropies = new DenseVector[Double](dim)
        (0 until dim).map(i => entropies(i) = entropyEstimator(X(i,::).t))


        // at most 10 sweeps with Jacoby rotations J(i,j,alpha) over all coordinate pairs (i,j) with i<j
        // target: maximum entropy = minimum contrast function
        var sweep = 0
        var contrastFcn = sum(entropies)
        if(verbose) System.out.print("\nContrast function at start: "+MathTools.round(contrastFcn,2)+"\n")
        var contrastFcn_new = contrastFcn
        var delta_CF = Double.MinValue                                  // to enter the loop
        while(sweep<10 & delta_CF < -Math.abs(0.01*contrastFcn)){

            var i=0
            while(i<dim){
                var j=i+1
                val row_i:DenseVector[Double] = X(i,::).t        // row_i(X)
                while(j<dim){

                    val row_j:DenseVector[Double] = X(j,::).t      // row_j(X)

                    // parallelize on angle search
                    val nThread: Int = Runtime.getRuntime().availableProcessors()
                    val futures = ListBuffer.empty[Future[(Double,Double)]]
                    for(p <- 0 until nThread){

                        futures += parallelAngleSearch(i,j,dPhi,p,nThread,X,entropies)
                    }
                    var dCF_opt = Double.MaxValue    // change in contrast function, wanted: large, negative
                    var phi_opt = 0.0
                    for(fut <- futures){

                        val (phi:Double,dCF:Double) = Await.result(fut, 10000 second)
                        if(dCF < dCF_opt){ dCF_opt = dCF; phi_opt=phi}
                    }
                    if(dCF_opt<0) {
                        // optimal angle and entropy decrease for coordinates (i,j) are set,
                        // apply rotation and update state to current optimum
                        val a = Math.cos(phi_opt)
                        val b = Math.sin(phi_opt)
                        val u_i = a*row_i-b*row_j
                        val u_j = b*row_i+a*row_j
                        X(i,::) := u_i.t
                        X(j,::) := u_j.t
                        entropies(i) = entropyEstimator(u_i)
                        entropies(j) = entropyEstimator(u_j)
                        rot.addRotation(i,j,phi_opt)
                    }
                    j+=1
                } // while j
                i+=1
            } // while i
            contrastFcn_new = sum(entropies)
            delta_CF = contrastFcn_new-contrastFcn
            if(vrbose){
                val distDecreasePcnt = MathTools.round(100*delta_CF/contrastFcn,2)
                var msg:String = "Sweep "+sweep+": contrast function = "+MathTools.round(contrastFcn_new,2)
                msg += ", decrease (%) = "+distDecreasePcnt
                System.out.println(msg)
            }
            contrastFcn = contrastFcn_new
            sweep+=1
        } // while sweep
        rot
    } // findOptimalRotationPar

    /** Search for counterclockwise optimal rotation phi in coordinate plane in coordinate i,j
    * (spanned by unit vectors e_i,e_j) starting at angle phi0 and stepping with step width dPhi up
        * to pi/4.
    *
    * @param X workspace for rotation of whitened data.
    * @param entropies: vector of entropies of the coordinates (rows) of X.
    * @return Ordered pair (phi_opt,dCF) where phi is the optimal angle of rotation and dCF the (additive)
    *         change in contrast function when ths rotation is applied to the data.
        */
    def optimalRotation(
        i:Int, j:Int, phi0:Double, dPhi:Double, X:DenseMatrix[Double], entropies:DenseVector[Double]
    ):(Double,Double) = {

        var phi = phi0
        var phi_opt=phi
        var dCF = 0.0                                // change in contrast function
        var dCF_opt = Double.MaxValue                // optimal change in contrast function
        val row_i:DenseVector[Double] = X(i,::).t      // row_i(X)
        val row_j:DenseVector[Double] = X(j,::).t      // row_j(X)

        while(phi<pi/4+1e-6){

            val a = Math.cos(phi)
            val b = Math.sin(phi)
            // try rotation J(i,j,phi)
            val u_i = a*row_i-b*row_j
            val u_j = b*row_i+a*row_j
            val H_i = entropyEstimator(u_i)
            val H_j = entropyEstimator(u_j)
            // change in contrast function if we apply this rotation
            // (X(i,::) -> u_i, X(j,::) -> u_j):
            val dCF = (H_i+H_j)-(entropies(i)+entropies(j))
            if(dCF < dCF_opt){ phi_opt = phi; dCF_opt = dCF }
            phi+=dPhi
        }
        (phi_opt,dCF_opt)
    }

    /** Search for angle phi_opt which minimizes the change dCF in contrast function (maximal decrease)
      * when rotating the data with the Jacobi rotation J(i,j,phi_opt). This rotates counterclockwise in the
      * coordinate plane in coordinates i,j (spanned by the standard unit vectors e_i,e_j).
      *
      * The overall search grid is equidistant in the interval [-pi/4,pi/4] with step width dPhi = pi/(2*K),
      * i.e. K equidistant points. The future number p searches all angles
      *
      *    phi = -pi/4 + j*dPhi, where j mod nThreads = p.
      *
      * @param i  coordinate of Jacobi rotation J(i,j,phi)
      * @param j  coordinate of Jacobi rotation J(i,j,phi)
      * @param p  angle search starts at -pi/4+p*phi_step
      * @param nThread number of threads employed.
      * @param dPhi  step width of search grid in [-pi/4,pi/4].
      * @param X data (row_i is the time series of coordinate X_i).
      * @param entropies: entropies(i) is the entropy H_i of X_i.
      * @return   Future with value (phi_opt,dCF), where dCF is the change in contrast function
      *           when rotating the data with the Jacobi rotation J(i,j,phi_opt)
      */
    def parallelAngleSearch(
        i:Int,j:Int,dPhi:Double,p:Int,nThread:Int,X:DenseMatrix[Double], entropies:DenseVector[Double]
    ): Future[(Double,Double)] = Future {

        val phi0 = -pi/4+p*dPhi
        optimalRotation(i,j,phi0,nThread*dPhi,X,entropies)
    }
}


object RadicalICA {

    val rnorm = new breeze.stats.distributions.Gaussian(0,1)
    /**
      * Factory function for allocating RadicalICA objects.
      *
      * @param data sample of multidimensional distribution, each column is one random value of the
      *              distribution.
      * @param nAngles number of equidistant grid points in [-pi/4,pi/4] scanned for optimal angles alpha
      *                of the Jacobi rotations.
      * @param entropyEstimator estimator for the marginal entropies. Provided are
      *      [[MathTools.entropyEstimatorVasicek]], [[MathTools.entropyEstimatorEmpiricalAdapted]] and
      *      [[MathTools.entropyEstimatorEmpiricalFast]]
      * listed from slowest and most precise to fastest and least precise. The last estimator can fail badly on
      * strongly clustered data (spiky densities).
      * @param doWhitenData padded data will be whitened before proceeding. Should always be done.
      *      But for testing purposes we want to be able to switch this off
      * @param doParallelSearch use the multithreaded version of the rotation search,
      *      the sequential version exists only for testing and debugging.
      */
    def apply(   data:DenseMatrix[Double],
                 nAngles:Integer,
                 entropyEstimator: (DenseVector[Double])=>Double,
                 doWhitenData:Boolean,
                 doParallelSearch:Boolean,
                 verbose:Boolean
    ) = new RadicalICA(data, nAngles, entropyEstimator, doWhitenData, doParallelSearch, verbose)

    def whiten(data:DenseMatrix[Double]):DenseMatrix[Double] = MatrixFunction.whiten(data)

    /** Repeat each column of data nPad times, each time adding an independent N(0,sigma²)
      * perturbation to each coordinate.
      * In other words the data sample distribution is convolved with multinormal N(0,sigma²I),
      * where I is the identity matrix in the same dimension as the data.
      *
      * @return padded data matrix consisting of the perturbed columns.
      */
    def pad(data:DenseMatrix[Double],nPad:Integer,sigma:Double):DenseMatrix[Double] = {

        val nSample = data.cols
        val dim = data.rows
        val paddedData = DenseMatrix.zeros[Double](dim,nPad*nSample)

        var row=0
        while(row<dim){
            var col=0
            while(col<nPad*nSample){
                paddedData(row,col)=data(row,col%nSample)+sigma*rnorm.draw()
                col+=1
            }
            row+=1
        }
        paddedData
    }

}
