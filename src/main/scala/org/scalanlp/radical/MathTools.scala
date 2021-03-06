package org.scalanlp.radical


import breeze.linalg.{max,DenseMatrix, DenseVector}
import breeze.numerics._


/**
 * Created by MeyerM on 03.12.2014.
 */
object MathTools {


    val pi = 3.1415926535897932

    /** Beta(a,b)=gamma(a)gamma(b)/gamma(a+b).*/
    def betaFcn(a:Double,b:Double) = Math.exp(lgamma(a))*Math.exp(lgamma(b))/Math.exp(lgamma(a+b))

    /** Vasicek entropy estimator for one dimensional random variable, see docs/RADICAL_ICA.pdf, p3, eq(12)).
      * This the slowest, most accurate estimator.
      *
      * @param sample: sample of random variable.
      * @param m spacing of Vasicek entropy estimator, see docs/RADICAL_ICA.pdf, p3, eq(12).
      */
    def entropyEstimateVasicek(sample:DenseVector[Double],m:Integer):Double = {

        val n=sample.length
        assert(m<n)
        val z=sample.toScalaVector().sortWith(_<_)     // ordered sample
        var sum=0.0; var i=0
        while(i+m<n){ sum+=Math.log(z(i+m)-z(i)); i+=1 }
        sum/(n-m)+Math.log((n+1).toDouble/m)
    }
    /** Convenient type modification of [[entropyEstimateVasicek]].
      * The slowest most accurate estimator on clustered data.
      */
    def entropyEstimatorVasicek(m:Int): (DenseVector[Double])=>Double =
        (sample:DenseVector[Double]) => entropyEstimateVasicek(sample,m)

    /** Entropy computed from histogram with 1+sqrt(sample size) bins of equal width.
      * This is very fast but may be very inaccurate on strongly clustered data (spiky densities),
      * see MathTests.testEntropySpikyDist.
      *
      * @param sample sample of random variable.
      */
    def entropyEstimateEmpiricalFast(sample:DenseVector[Double]):Double = {

        val hist = NaiveHistogram(sample)
        hist.entropy
    }
    /** Convenient type modification of [[entropyEstimateEmpiricalFast]].*/
    def entropyEstimatorEmpiricalFast:(DenseVector[Double])=>Double =
        (sample:DenseVector[Double]) => entropyEstimateEmpiricalFast(sample)
    /** Entropy computed from adapted histogram with nBins bins (approximately).
      * Much better (and slower) than [[entropyEstimatorEmpiricalFast]] on clustered data and faster,
      * but not as good as [[entropyEstimatorVasicek]].
      *
      * @param sample sample of random variable.
      * @param nBins number of bins in the histogram used for estimation.
      */
    def entropyEstimateEmpiricalAdapted(sample:DenseVector[Double],nBins:Int):Double = {

        val hist = AdaptedHistogram(sample,nBins)
        hist.entropy
    }
    /** Convenient type modification of [[entropyEstimateEmpiricalFast]].*/
    def entropyEstimatorEmpiricalAdapted(nBins:Int):(DenseVector[Double])=>Double =
        (sample:DenseVector[Double]) => entropyEstimateEmpiricalAdapted(sample,nBins)

    /**
      * Entropy of the uniform distribution on [0,a].
      */
    def entropyUniform(a:Double) = Math.log(a)
    /**
      * Variance of the uniform distribution on [0,a].
      */
    def varianceUniform(a:Double) = a*a/12
    /**
      * Mean of the uniform distribution on [0,a].
      */
    def meanUniform(a:Double) = a/2



    /**
      * Entropy of the Gaussian distribution N(mu,sigma).
      * Note that parametrization is _not_ N(mu,sigma*sigma) as usual.
      *
      * @param mu: mean
      * @param sigma: standard deviation
      */
    def entropyGaussian(mu:Double, sigma:Double) = (1+Math.log(2*pi*sigma*sigma))/2



    /**
      * Entropy of the Exponential distribution \[E(\lambda)\] (mean is \[1/\lambda\]).
      *
      * @param lambda: rate
      */
    def entropyExponential(lambda:Double) = 1-Math.log(lambda)
    /**
      * Variance of the Exponential distribution \[E(\lambda)\] (mean is \[1/\lambda\]).
      *
      * @param lambda: rate
      */
    def varianceExponential(lambda:Double) = 1.0/(lambda*lambda)
    /**
      * Mean of the Exponential distribution \[E(\lambda)\] (mean is \[1/\lambda\]).
      *
      * @param lambda: rate
      */
    def meanExponential(lambda:Double) = 1.0/lambda



    /**
      * Entropy of the Beta distribution Beta(a,b) on [0,1].
      *
      * @param a: parameter usually called alpha
      * @param b: parameter usually called beta
      */
    def entropyBeta(a:Double,b:Double) =
        lgamma(a)+lgamma(b)-lgamma(a+b) - (a-1)*digamma(a) - (b-1)*digamma(b) + (a+b-2)*digamma(a+b)
    /**
      * Variance of the Beta distribution Beta(a,b) on [0,1].
      *
      * @param a: parameter usually called alpha
      * @param b: parameter usually called beta
      */
    def varianceBeta(a:Double,b:Double) = a*b/((a+b)*(a+b)*(a+b+1))
    /**
      * Mean of the Beta distribution Beta(a,b) on [0,1].
      *
      * @param a: parameter usually called alpha
      * @param b: parameter usually called beta
      */
    def meanBeta(a:Double,b:Double) = a/(a+b)




    /**
      * Entropy of the Gamma distribution Gamma(a,b) on [0,+oo).
      *
      * @param shape: shape parameter
      * @param scale: scale parameter
      */
    def entropyGamma(shape:Double,scale:Double) = shape + log(scale) + lgamma(shape) + (1-shape)*digamma(shape)
    /**
      * Variance of the Gamma distribution Gamma(a,b) on [0,+oo).
      *
      * @param shape: shape parameter
      * @param scale: scale parameter
      */
    def varianceGamma(shape:Double,scale:Double) = shape*scale*scale
    /**
      * Mean of the Gamma distribution Gamma(a,b) on [0,+oo).
      *
      * @param shape: shape parameter
      * @param scale: scale parameter
      */
    def meanGamma(shape:Double,scale:Double) = shape*scale


    /** @return value y=f(x) from linear interpolation of
      *  values (x1,y1=f(x1)), (x2,y2=f(x2))
      */
    @inline
    final def linInt(x:Double,x1:Double,y1:Double,x2:Double,y2:Double) =
        if      (abs(x-x1)<1e-10) y1
        else if (abs(x-x2)<1e-10) y2
        else ((x-x1)*y1+(x2-x)*y2)/(x2-x1)

    /** Sum x(0)+x(1)+...+x(upToIdx) */
    @inline
    final def sumFromLeft(x:Seq[Double], upToIdx: Int) = {

        var sum = 0.0; var i=0
        while(i<upToIdx){ sum+=x(i); i+=1 }
        sum
    }
    /** Apply function f to each entry of M, return new matrix.*/
    final def applyFcn(u:DenseVector[Double],f:Double=>Double):DenseVector[Double] = {

        val res = new DenseVector[Double](u.length)
        var i=0
        while(i<u.length){ res(i)=f(u(i)); i+=1 }
        res
    }
    /** apply function f to each entry of M, return new matrix. */
    final def applyFcn(M:DenseMatrix[Double],f:Double=>Double):DenseMatrix[Double] = {

        val res = new DenseMatrix[Double](M.rows,M.cols)
        var row=0
        while(row<M.rows){
            var col=0
            while(col<M.cols){ res(row,col)=f(M(row,col)); col+=1 }
            row+=1
        }
        res
    }
    final def round(u:Double,decimals:Int):Double = {
        val f=(1 until decimals).foldLeft(10.0)((v:Double,k:Int)=>10*v)   // 10^decimals
        Math.floor(u*f)/f
    }
    final def round(v:DenseVector[Double],decimals:Int):DenseVector[Double] = applyFcn(v,(u:Double)=>round(u,decimals))
    final def round(M:DenseMatrix[Double],decimals:Int):DenseMatrix[Double] = applyFcn(M,(u:Double)=>round(u,decimals))

    /** The fastest way to get sample minimum and maximum (3n/2 comparisons).
      *
      * @param sample data sample, each column is one realization of the random vector.
      * @return pair (min(sample),max(sample)).
      */
    final def extremeValues(sample:DenseVector[Double]):(Double,Double) = {

        val n = sample.length
        var min = Double.MaxValue
        var max = Double.MinValue
        var i=0
        while(2*i+1<n){

            val a = sample(2*i)
            val b = sample(2*i+1)
            if(a < b) {
                if (a < min) min = a
                if (b > max) max = b
            } else { // b<=a
                if (b < min) min = b
                if (a > max) max = a
            }
            i+=1
        }
        // if sample length is even we have not checked the last element
        if(2*i==n){
            val a = sample(n-1)
            if(a < min) min=a else if (a>max) max=a
        }
        (min,max)
    }
    /** Check if U'U yields the identity matrix I.
      * @return max(|U'U-I|) < tolerance.*/
    def isOrthogonal(U:DenseMatrix[Double],tolerance:Double): Boolean = {

        assert(U.rows==U.cols,"Matrix U is not square.\n")
        val I = DenseMatrix.eye[Double](U.rows)
        val C = U.t * U
        val diff = max(abs(C-I))
        diff < tolerance
    }
}
