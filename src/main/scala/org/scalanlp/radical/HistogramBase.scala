package org.scalanlp.radical

import net.jafama.FastMath

/**
  * Created by oar on 7/2/16.
  * This contains all the methods which can be implemented on the basis of the methods in the
  * interface Histogram.
  * The specialized histogram classes then contain indiosyncratic methods.
  */
abstract class HistogramBase extends Histogram {

    /**
      * This is the integral \[\int f(x)log(f(x))dx\] where the density f(x) is approximated by the stepfunction
      * implied by the bin counts.
      *
      * @return
      */
    def entropy:Double = if (numBins==1) 0.0 else {

        var sum=0.0
        var i = 0
        while(i < numBins){
            val b = bin(i)
            val binWidth = b._2-b._1
            // density = normalized bin height
            val f_i = binCount(i)/(sampleSize*binWidth)
            if(f_i>1e-12) sum -= f_i*FastMath.log(f_i)*binWidth
            i+=1
        }
        sum
    }
}
