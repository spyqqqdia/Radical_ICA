package org.scalanlp.radical

import breeze.linalg.{Axis, DenseMatrix, DenseVector, eigSym, sum, svd, max => bn_max}
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
  * Vasicek estimator. The class is adapted to the use in the ICA algorithm. This explains the data structures
  * chosen.
  */
class Histogram(val sample:DenseVector[Double]) {

    val extremes = MathTools.extremeValues(sample)
    val min = extremes._1
    val max = extremes._2
    // if max-min < 1e-12 we'll view the random variable as constant
    val nBins:Int = if(max-min<1e-12) 1 else (1+FastMath.sqrt(sample.length)).toInt   // number of bins
    val h = (max-min)/nBins     // bin width
    // compute bin counts
    val bins = DenseVector.zeros[Double](nBins)
    if(nBins==1) bins(0)=sample.length else
    {
        val d = 1.0/h
        var i=0;
        while(i<sample.length){

            var j = ((sample(i)-min)*d).toInt
            if(j>=nBins) j=nBins-1       // when the max is hit
            bins(j)+=1;
            i+=1;
        }
    }

    /**
      * This is the integral \[\int f(x)log(f(x))dx\] where the density f(x) is approximated by the stepfunction
      * implied by the bin counts.
      * @return
      */
    def entropy:Double = if(nBins==1) 0.0 else {

        val d = 1.0/h
        var sum=0.0
        var i=0
        while(i<nBins){
            // density = normalized bin height (1 = area = h*(sum of bin heights))
            val f_i = d*bins(i)/sample.length
            if(f_i>1e-12) sum -= f_i*FastMath.log(f_i)
            i+=1
        }
        h*sum
    }
}



object Histogram {

    /** Factory method.*/
    def apply(sample:DenseVector[Double]) = new Histogram(sample)
}
