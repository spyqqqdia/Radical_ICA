package org.scalanlp.radical

/**
  * Created by oar on 7/2/16.
  * Class to represent the empirical distribution of a sample. For this we need only the probabilities
  * p_j=P(X=x_j) and bin boundaries but not the values x_j themselves which therefore will be disregarded.
  *
  * The class HeavyHistogram implements a histogram that actually allocates bins and puts the sample values
  * in the respective bins. This is useful for bin sorts.
  */
trait Histogram {

    /** Pair (xl,xr) where the i-th bin contains all values x of the sample with
      * xl <= x < xr except for the last bin where the inequality is xl <= x <= xr.
      *
      * @param i number of bin, count starts at zero, i=0,...,numBins-1.
      * */
    def bin(i:Int):(Double,Double)
    /** Number of bins. */
    def numBins:Int
    /** Number of elements in i-th bin, i=0,1,...,numBins-1 */
    def binCount(i:Int):Int
    def sampleSize:Int
}
