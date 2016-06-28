and forces us to Title:  Linux Administration
CSS:    css/bootstrap_min.css
CSS:    css/tiny.css
Author: Michael J. Meyer  
Date:   June 24, 2016



#Radical ICA#

This project implements the independent component analysis (`ICA`) algorithm using Scala 2.11.8 and JRE 1.8.0_76 with IntelliJ 2016.1.2.   
Original authors of the algorithm: Eric G. Learned-Miller, John W. Fisher III.
 
* [Introduction](radical_ica.pdf)
* [Radical ICA, original paper](RadicalICA_MilFish.pdf)
* [Scaladoc](index.html)
* [Prof. Learned-Miller, homepage](https://people.cs.umass.edu/~elm/)
* [Prof. Learned-Miller, R, Matlab code](https://people.cs.umass.edu/~elm/papers_by_code.html)
* [Radical ICA hompage at umass](https://people.cs.umass.edu/~elm/ICA/)

Independent component analysis is a search for a rotation of whitened multidimensional data (covariance matrix is the identity matrix)  
which increases the independence of the components (coordinates). 

The rotation is chosen so as to minimize a contrast function. The contrast function is the sum of the marginal entropies  
(entropies of the components) of the sample plus a constant. Therefore an entropy estimator has to be employed.  
The original algorithm uses the Vasicek entropy estimator based on the order statistic which makes it necessary to sort the components  
of the sample. This is O(n*log(n)) in the sample size n.

The Radical ICA algorithm performs an exhaustive search and is therefore naturally slow. Here an attempt is made to speed this up  
as much as possible. To this end a parallelized implementation is provided. In addition the empirical entropy estimator is provided  
as an alternative to the Vasicek estimator. The algorithm is

* O(dim²nlog(n)) with the Vasicek entropy estimator and
* O(dim²n) with the empirical entropy estimator,

where dim is the dimension of the data. With a sample of size 30000 in dimension 6 the version using the empirical entropy estimator  
is more than 10 times faster (processor: core i5 5765C with 128MB L3 cache).

The two entropy estimators have been tested on a variety of distributions as follows: 1000 samples of size 30000 and 3000 are  
generated and the smaller sample padded to size 30000 as described in the original paper. The entropies of the samples are  
then computed with both estimators and mean and variance of the estimate computed. These values are then compared with the  
true entropy of the distribution generating the samples. In all cases tested the two estimators show similar performance,  
see [results](../results/EntropyEstimation.csv). 


