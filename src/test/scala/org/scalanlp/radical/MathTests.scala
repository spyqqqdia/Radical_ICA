package org.scalanlp.radical

import breeze.linalg.DenseVector
import breeze.stats.{mean, variance}
import breeze.stats.distributions._


/**
  * Created by oar on 6/22/16.
  */
object MathTests {

    // see function entropyTestForDistribution
    type entropyTestResult = (                                        // pair of tuples
                    (String,String,Double,Double,Double,Double,Double),
                    (String,String,Double,Double,Double,Double,Double)    )

    /** Helper function to format the output of function [[MathTests.entropyTestForDistribution]] below.
      *
      * @param dist distribution
      * @param result return value of [[MathTests.entropyTestForDistribution]]
      * @return formatted string for export as csv file.
      */
    def formatEntropyTestResult(dist:ContinuousDistr[Double], result:entropyTestResult):String = {

        val rd = (u:Double) => MathTools.round(u,4)
        var str:String = ""
        val rV = result._1          // results for Vasicek estimator
        str += rV._1+";  "+rV._2+";  "+rd(rV._3)+";  "+rd(rV._4)+";  "+rd(rV._5)+";  "+rd(rV._6)+";  "+rd(rV._7)+"\n"
        val rE = result._2          // results for empirical estimator
        str += rE._1+";  "+rE._2+";  "+rd(rE._3)+";  "+rd(rE._4)+";  "+rd(rE._5)+";  "+rd(rE._6)+";  "+rd(rE._7)+"\n"
        str
    }


    /**
      * Generates 1000 samples of of size 30000 and 1000 samples of size 3000 padded with multiplicity 10
      * (so that the size of the padded sample (i.e. sample convolved with centered Gaussian) is again 30000).
      * Then computes the Vasicek and empirical entropy estimates of both samples and reports mean and standard
      * deviation of the estimates.
      *
      * The computation is carried out separately for the Vasicek and empirical estimates so that we can time the
      * estimators. In our examples the empirical estimator performs about as well as the Vasicek estimator but is
      * 10 times faster. This is masked partially by the sample generation so that only a speedup of 4-5 is observed.
      *
      * The reason for this is that the Vasicek estimator uses the order statistic which needs a sort of O(n*log(n))
      * whereas the empirical estimator is O(n) in the sample size n.
      *
      * @param dist: a continuous distribution
      * @param entropy: entropy of the distribution
      * @param sigmaPad: standard deviation of normal distribution used in padding (convolution),
      *                  see [[DataGenerator.smoothOneDimensionalData]]
      * @return pair (resultsVasicek,resultsEmpirical) where
      *     resultsVasicek = (
      *        cDist.toString(), "Vasicek", entropy(cDist), meanEstimate, sdEstimate, meanEstimatePadded, sdEstimatePadded
      *     )
      *     resultsEmpirical is analogous (here "sd" means standard deviation).
      */
    def entropyTestForDistribution(dist:ContinuousDistr[Double],entropy:Double,sigmaPad:Double):entropyTestResult = {

        var msg  = "\nDistribution: "+dist
        System.out.println(msg)

        val sampleSize = 30000
        val nTests = 1000
        val entropies = DenseVector.zeros[Double](nTests)
        val entropiesP = DenseVector.zeros[Double](nTests)   // padded sample

        // Vasicek
        Timer.start
        (0 until nTests).map(i => {

            val s = dist.sample(sampleSize)
            val sample = DenseVector.tabulate(sampleSize){i => s(i)}    // turn s into vector
            val samplePadded = DataGenerator.smoothOneDimensionalData(sample(0 until (sampleSize/10)),10,sigmaPad)
            val m = Math.sqrt(sampleSize).toInt
            entropies(i) = MathTools.entropyEstimateVasicek(sample,m)
            entropiesP(i) = MathTools.entropyEstimateVasicek(samplePadded,m)

        })
        Timer.stop
        msg = "Vasicek entropy, "+Timer.report
        System.out.println(msg)

        var mu = mean(entropies);
        mu = MathTools.round(mu,4)
        var stdv = Math.sqrt(variance(entropies))
        stdv = MathTools.round(stdv,4)
        var muP = mean(entropiesP)
        muP = MathTools.round(muP,4)
        var stdvP = Math.sqrt(variance(entropiesP))
        stdvP = MathTools.round(stdvP,4)

        val resVasicek = (dist.toString(), "Vasicek", entropy,mu,stdv,muP,stdvP)

        // empirical
        Timer.start
        (0 until nTests).map(i => {

            val s = dist.sample(sampleSize)
            val sample = DenseVector.tabulate(sampleSize){i => s(i)}    // turn s into vector
            val samplePadded = DataGenerator.smoothOneDimensionalData(sample(0 until (sampleSize/10)),10,sigmaPad)
            entropies(i) = MathTools.entropyEstimateEmpiricalFast(sample)
            entropiesP(i) = MathTools.entropyEstimateEmpiricalFast(samplePadded)
        })
        Timer.stop
        msg = "Empirical entropy, "+Timer.report
        System.out.println(msg)

        mu = mean(entropies);
        mu = MathTools.round(mu,4)
        stdv = Math.sqrt(variance(entropies))
        stdv = MathTools.round(stdv,4)
        muP = mean(entropiesP)
        muP = MathTools.round(muP,4)
        stdvP = Math.sqrt(variance(entropiesP))
        stdvP = MathTools.round(stdvP,4)

        val resEmpirical = (dist.toString(), "Empirical", entropy,mu,stdv,muP,stdvP)

        (resVasicek,resEmpirical)
    }
    /** Carries out the [[entropyTestForDistribution]] of the entropy estimators for
      * Uniform(0,10), Gaussian(0,5) and Exponential(10) as well as Gamma(2,2) and Gamma(7,1) and
      * Beta(0.5,0.5), Beta(1,3), Beta(5,1), and Beta(2,5)  distributions. See en.wikipedia for the shape
      * of the Gamma and Beta distributions.
      */
    def entropyTest():Unit = {

        val doUniform = true
        val doGaussian = true
        val doExponential = true
        val doGamma = true
        val doBeta = true

        var msg = "Vasicek versus empirical entropy, padding sample with sigma = SD(distribution)/10."
        System.out.println(msg)
        var result:entropyTestResult = null
        // csv header line
        var results:String = "Distribution;  Estimator;  Entropy;  Mean_estimate;  SD_estimate;  "
        results += "Mean_estimate_padded;  SD_estimate_padded\n"


        // Uniform distribution on [0.a]
        var a = 10.0
        var sigmaPad = Math.sqrt(MathTools.varianceUniform(a))/10      // SD/7 seems to work
        var cDist:ContinuousDistr[Double] = new Uniform(0.0,a)
        var entropy = MathTools.entropyUniform(a)
        if(doUniform){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }


        // Normal distribution
        val mu = 0.0
        val sigma = 5.0
        sigmaPad = sigma/7
        cDist = new Gaussian(mu,sigma)
        entropy = MathTools.entropyGaussian(mu,sigma)
        if(doGaussian){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }


        // Exponential dist, sigmaPad must be set very small!
        val lambda = 10.0
        sigmaPad = Math.sqrt(MathTools.varianceExponential(lambda))/10
        cDist = new Exponential(lambda)
        entropy = MathTools.entropyExponential(lambda)
        if(doExponential){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }


        // Gamma dist, two dissimilar shapes, see en.wikipedia
        var shape=2.0; var scale=2.0
        sigmaPad = Math.sqrt(MathTools.varianceGamma(shape,scale))/10
        cDist = new Gamma(shape,scale)
        entropy = MathTools.entropyGamma(shape,scale)
        if(doGamma){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }

        shape=7.0; scale=1.0
        sigmaPad = Math.sqrt(MathTools.varianceGamma(shape,scale))/10
        cDist = new Gamma(shape,scale)
        entropy = MathTools.entropyGamma(shape,scale)
        if(doGamma){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }


        // Beta dist, four dissimilar shapes, see en.wikipedia
        a=0.5; var b=0.5;
        cDist = new Beta(a,b)
        sigmaPad = Math.sqrt(MathTools.varianceBeta(a,b))/10
        entropy = MathTools.entropyBeta(a,b)
        if(doBeta){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }

        a=1.0; b=3.0;
        cDist = new Beta(a,b)
        sigmaPad = Math.sqrt(MathTools.varianceBeta(a,b))/10
        entropy = MathTools.entropyBeta(a,b)
        if(doBeta){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }

        a=5.0; b=1.0;
        cDist = new Beta(a,b)
        sigmaPad = Math.sqrt(MathTools.varianceBeta(a,b))/10
        entropy = MathTools.entropyBeta(a,b)
        if(doBeta){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }

        a=2.0; b=5.0;
        cDist = new Beta(a,b)
        sigmaPad = Math.sqrt(MathTools.varianceBeta(a,b))/10
        entropy = MathTools.entropyBeta(a,b)
        if(doBeta){

            result = entropyTestForDistribution(cDist, entropy, sigmaPad)
            results += formatEntropyTestResult(cDist,result)
        }
        System.out.println(results)
        System.out.println("\nFinished.")
    }

    /** Comparing the entropy estimators on a sample of size sampleSize for the following distribution P:
      * P = 0.5*Unif(0,1)+0.5*Unif(0,0.0001) with density
      * f(u) = 0.5*(1+1/0.0001) on [0,0.0001] and f(u)=0.5 on the remaining interval.
      * The spike at the left endpoint is a problem for the naive estimator with equal width bins.
      */
    def testEntropySpikyDist(sampleSize:Int): Unit = {

        val distributions = Vector[ContinuousDistr[Double]](Uniform(0,1d), Uniform(0,0.0001d))
        val weights = Vector[Double](0.5, 0.5)

        val P:ContinuousDistr[Double] = ContinuousMixtureDistribution(distributions,weights)
        val sample: DenseVector[Double] = new DenseVector[Double](P.sample(sampleSize).toArray)

        val entHistFast = MathTools.entropyEstimateEmpiricalFast(sample)
        val entHistAdapted = MathTools.entropyEstimateEmpiricalAdapted(sample,sample.size/50)
        val m = Math.sqrt(sampleSize).toInt    // spacing
        val entVasicek = MathTools.entropyEstimateVasicek(sample,m)
        val entropy = -(0.0001*5000.5*Math.log(5000.5)+0.9999*0.5*Math.log(0.5))

        var message = "Entropy, theoretical: "+MathTools.round(entropy,4)+",\n"
        message += "entropy, naive histogram estimator: "+MathTools.round(entHistFast,4)+",\n"
        message += "entropy, adapted histogram estimator estimator: "+MathTools.round(entHistAdapted,4)+",\n"
        message += "entropy, Vasicek estimator: "+MathTools.round(entVasicek,4)+".\n"

        print(message)

    }


}
