package org.scalanlp.radical

import breeze.linalg.{DenseMatrix, DenseVector}
import breeze.stats.distributions._


/**
  * Created by oar on 7/3/16.
  *
  * Generates samples of centered multidimensional data with independent components rotated by a given rotation
  * and padded to larger size (and smoothed) by convolution with independent Gaussian noise.
  *
  * These are the test samples we will use for our Radical ICA algorithm.
  */
object DataGenerator {

    /** N(0,1) generator.*/
    val rnorm = new breeze.stats.distributions.Gaussian(0,1)

    /** Standard list of distributions for the coordinates. Some of these are moderately somewhat clustered:
      *
      **/
    val distributionList = List[ContinuousDistr[Double]](

        Uniform(0,1), Exponential(5), Gaussian(0,0.25), Gamma(2,2), Gamma(1,2), BiExponential(10)
    )

    /** Repeat each entry in sample nPad ties while adding a N(0,sigma²) draw to it.
      *
      * This means that we generate a sample of size nPad*sample.size of the convolution of the distribution
      * generating the sample with the N(0,sigma²) distribution (Gaussian smoothing).
      *
      */
    def smoothOneDimensionalData(sample:DenseVector[Double], nPad:Int, sigma:Double):DenseVector[Double] = {

        val n = sample.length*nPad
        val res = DenseVector.zeros[Double](n)
        var i=0
        for(i <- 0 until n) res(i) = sample(i%sample.length)+sigma*rnorm.draw()
        res
    }
    /** Each column of sample is interpreted as a realization of a random vector. The function enlarges the
      * sample by repeating each column nPad times and adding N(0,sigma²) variates to each coordinate.
      *
      * This means is that we we generate a sample of size nPad*sample.size of the convolution of the distribution
      * generating the sample with the Gaussian distribution N(0,sigma²I) in the same dimension as the sample data
      * (Gaussian smoothing).
      *
      */
    def smoothMultiDimensionalData(sample:DenseMatrix[Double], nPad:Int, sigma:Double):DenseMatrix[Double] = {

        val rows = sample.rows
        val cols = sample.cols
        val colsPadded = sample.cols*nPad
        val res = DenseMatrix.zeros[Double](rows,colsPadded)
        for(i <- 0 until rows; j <- 0 until colsPadded) res(i,j) = sample(i,j%cols)+sigma*rnorm.draw()
        res
    }


    /** Generates a data sample Z in dimension dim = coordinateDistributions.size as follows:
      * (a) Let X be a random vector in dimension dim with independent coordinates X_i with distribution
      *     coordinateDistributions(i) each.
      * (b) A sample of size sampleSize/40 is drawn from X and padded to size sampleSize (padding factor 40)
      *     using the function [[smoothMultiDimensionalData]] (smoothing by convolution with Gaussian noise).
      * (b) The coordinates of the padded sample are rescaled so that they have standard deviation 1. This yields a
      *     sample Y which is white (cov(Y)=I, the identity matrix).
      * (c) The sample Y is rotated by application of the rotation Q.
      *
      * @param sigma standard deviation of smoothing noise, useful default: 0.175.
      * @param rot rotation applied to the whitened data Y.
      * @return Z=QY
      */
    def rotatedWhitenedSample(
             coordinateDistributions:List[ContinuousDistr[Double]], sampleSize:Int, sigma:Double, rot:Rotation
    ):DenseMatrix[Double] = {

        val dim = coordinateDistributions.size
        if(dim!=rot.dim){

            val msg = "Dimension of rotation must equal the dimension "+dim+" of the data but is = "+rot.dim+".\n"
            throw new IllegalArgumentException(msg)
        }
        val unpaddedData = DenseMatrix.zeros[Double](dim,sampleSize/40)
        for(i <- 0 until dim; j <- 0 until sampleSize/40) unpaddedData(i,j)=coordinateDistributions(i).draw()

        val paddedData = smoothMultiDimensionalData(unpaddedData,40,sigma)
        rot()*paddedData
    }
}
