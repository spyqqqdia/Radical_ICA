package org.scalanlp.radical

import breeze.linalg.DenseVector
import net.jafama.FastMath

/**
  * Created by oar on 6/22/16.
  *
  * Class to represent the empirical distribution of a sample. With n denoting the sample size,
  * 1+sqrt(n) equal width bins will be allocated spanning the sample from minimum to maximum value.
  *
  * The intended applications is estimation of the sample entropy. For this we need only the probabilities
  * p_j=P(X=x_j) but not the values x_j (bin means) which therefore will be disregarded.
  *
  * The emphasis is on speed. In particular we circumvent the O(nlog(n)) sorting of the sample needed in the
  * Vasicek estimator. This will not work for highly clustered data (see MathTests.testEntropySpikyDist.
  * For such data the [[AdaptedHistogram]] is a better, but slower choice.
  */
class NaiveHistogram(val sample:DenseVector[Double]) extends HistogramBase {

    val extremes = MathTools.extremeValues(sample)
    val min = extremes._1
    val max = extremes._2
    // if max-min < 1e-12 we'll view the random variable as constant
    val numBins:Int = if(max-min<1e-12) 1 else (1+FastMath.sqrt(sample.length)).toInt   // number of bins
    val h = (max-min)/numBins     // bin width
    // compute bin counts
    val binCounts = DenseVector.zeros[Int](numBins)
    if(numBins==1) binCounts(0)=sample.length else
    {
        val d = 1.0/h
        var i=0
        while(i<sample.length){

            var j = ((sample(i)-min)*d).toInt
            if(j>=numBins) j = numBins-1       // when the max is hit
            binCounts(j)+=1;
            i+=1;
        }
    }
    def sampleSize = sample.size
    def binCount(i:Int) = binCounts(i)
    def bin(i:Int) = {

        assert((i>=0)&(i<=numBins), "Bin number "+i+" not in [0,"+(numBins-1)+"].\n")
        (i*h,(i+1)*h)
    }

}
object NaiveHistogram {

    /** Factory method.*/
    def apply(sample:DenseVector[Double]) = new NaiveHistogram(sample)
}
