package org.scalanlp.radical

import org.apache.commons.math3.random.MersenneTwister
import breeze.stats.distributions.{ContinuousDistr,RandBasis,Rand}

/**
  * Convex sum of continuous distributions \sum_jp_jDist_j in Euclidean space with weights p_j
  * (positive, summing to one).
  *
  * Generates draws as follows: Let N be the number of member distributions D_j. First you throw a loaded die with
  * N faces with probability p_j for face number j. If face number j comes up you draw a value from distribution
  * Dist_j. This is repeated for very single draw.
  *
  */
class ContinuousMixtureDistribution(val distributions:Vector[ContinuousDistr[Double]] , val weights:Vector[Double])
                                   (implicit basis: RandBasis = Rand)
extends ContinuousDistr[Double] {

    if(distributions.size != weights.size) {

        val message = "Number of distributions (" + distributions.size + ") not equal to number of weights (" + weights.size + ").\n"
        throw new IllegalArgumentException(message)
    }
    // should really use machine epsilon
    if(!(weights.sum==1.0)) {

        val message = "Weights not adding up tp one, weights.sum = "+weights.sum + ".\n"
        throw new IllegalArgumentException(message)
    }
    {
        val weightsArePos = weights.foldLeft(true)((c:Boolean,w:Double)=>(w>0.0)&c)
        if(!weightsArePos){

            val message = "Weights not all positive, weight:\n"+weights+ ".\n"
            throw new IllegalArgumentException(message)
        }
    }
    def numDists = weights.size
    // random number generator uniform in (0,1,...,numDists-1)
    val rng = basis.randInt(numDists)

    // density will already be defined in normalized form
    override def logNormalizer = 0.0
    /** This is the logarithm of the density (already normalized). */
    override def unnormalizedLogPdf(x:Double): Double = {

        // sorry, no scala idiomatic code since too slow
        var f = 0.0; var i=0
        while(i<numDists){

            val dist = distributions(i)
            f+= weights(i)*dist.unnormalizedPdf(x)/Math.exp(dist.logNormalizer)
            i+=1
        }
        Math.log(f)
    }
    override def draw():Double = {

        val i:Int = rng.get()
        distributions(i).draw
    }
}
object ContinuousMixtureDistribution {

    implicit val basis = new RandBasis(new MersenneTwister())
    /** Factory method.*/
    def apply(distributions:Vector[ContinuousDistr[Double]], weights:Vector[Double]) =
        new ContinuousMixtureDistribution(distributions,weights)
}
