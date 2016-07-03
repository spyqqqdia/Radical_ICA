package org.scalanlp.radical

import breeze.stats.distributions.{Exponential, ContinuousDistr, Rand, RandBasis}
import org.apache.commons.math3.random.MersenneTwister

/**
  * Created by oar on 7/3/16.
  *
  * Distribution on (-oo,+oo) with unnormalized density equal to the density of the Exponential distribution
  * plus this same density reflected about the origin. Samples as follows: draw a sample x from Exponential(rate),
  * then return x or -x with probability 0.5 each.
  *
  * @param rate rate parameter, must be >= 1e-15.
  *
  */
class BiExponential(rate:Double)(implicit basis: RandBasis = Rand) extends ContinuousDistr[Double] {

    val rng = basis.randInt(2) // unbiased coin
    val E = Exponential(rate)

    if(rate<1e-15){

        val msg = "Rate paramter must be >= 1e-15 but is = "+rate
        throw new IllegalArgumentException(msg)
    }

    // density will already be defined in normalized form
    override def logNormalizer = Math.log(2.0/rate)
    /** This is the logarithm of the density (already normalized). */
    override def unnormalizedLogPdf(x:Double): Double = -rate*Math.abs(x)
    override def draw():Double = {

        val x = E.draw()
        val i = rng.get()
        if(i==0) x else -x
    }
}
object BiExponential {

    implicit val basis = new RandBasis(new MersenneTwister())
    /** Factory method.*/
    def apply(rate:Double) = new BiExponential(rate)

}

