package org.scalanlp.radical

import scala.collection.immutable.VectorBuilder
import breeze.linalg.DenseVector
import breeze.stats.distributions.{Rand, RandBasis}
import net.jafama.FastMath
import org.apache.commons.math3.random.MersenneTwister

/**
  * Created by oar on 7/2/16.
  *
  * This is a histogram that adapts itself to the sample as follows: it draws a sub-sample of size nBins-1,
  * sorts the sub-sample and uses sample min, sample max and the sorted values of the sub-sample to set
  * the bin boundaries.
  *
  * The histogram may contain fewer bins than nBins, since duplicate values have to be removed from
  * the sub-sample.
  *
  * The intent is to have narrow bins in places where the data are clustered for better resolution.
  * For samples of size about 50000 setting nBins = sampleSize/50 can already handle some quite spiky
  * distributions, see MathTests.testEntropySpikyDist.
  *
  */
class AdaptedHistogram(val sample:DenseVector[Double], val nBins:Int)(implicit basis: RandBasis = Rand)
extends HistogramBase {

    val extremes = MathTools.extremeValues(sample)
    val min = extremes._1
    val max = extremes._2

    val rng = basis.uniform

    val binBoundaries:DoubleGrid = {

        val gridPoints = Vector.newBuilder[Double]
        gridPoints+=min; gridPoints+=max
        // fill with sub sample of size nBins-1
        var i=0
        while(i<nBins-1){

            val u = sampleSize*rng.draw()
            gridPoints+=sample(u.toInt)
            i+=1
        }
        val grid = gridPoints.result().sortWith(_<_).distinct
        DoubleGrid(grid)
    }


    // compute bin counts
    val binCounts = DenseVector.zeros[Int](numBins)
    if(numBins==1) binCounts(0)=sampleSize else
    {
        var i=0;
        while(i<sampleSize){

            var j = binBoundaries.leftIndex(sample(i))
            if(j== -1) j=0
            binCounts(j)+=1;
            i+=1;
        }
    }
    def numBins = binBoundaries.lastIdx
    def sampleSize = sample.size
    def binCount(i:Int) = binCounts(i)
    def bin(i:Int) = {

        assert((i>=0)&(i<=numBins), "Bin number "+i+" not in [0,"+(numBins-1)+"].\n")
        (binBoundaries(i),binBoundaries(i+1))
    }

}
object AdaptedHistogram {

    implicit val basis = new RandBasis(new MersenneTwister())
    /** Factory method.*/
    def apply(sample:DenseVector[Double],numBins:Int) = new AdaptedHistogram(sample,numBins)
}
