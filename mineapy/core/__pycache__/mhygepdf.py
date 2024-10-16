###  this function to calculate probability density function for multivariate hypergeometric function
import numpy as np
import scipy.special as sc
def mhygepdf(m,n):
    #MHYGEPDF Multivariate hypergeometric probability density function.
    # X = MHYGEPDF(M,N) returns the multivariate hypergeometric probability
    # density function at M with integer parameters in N. Note: The density
    # function is zero unless the elements in M are integers.
    #
    # If there are m_i elements of class i in a population and you take n
    # elements at random without replacement, then the number of elements of
    # each class in the sample (x_1,x_2,...,x_c) has the multivariate
    # hypergeometric distribution. This has the same relationship to the
    # multinomial distribution that the hypergeometric distribution has to the
    # binomial distribution--the multinomial distribution is the
    # 'with-replacement' distribution and the multivariate hypergeometric is
    # the 'without-replacement' distribution.
    #
    # Syntax: function x = mhygepdf(m,n)
    #
    # Inputs:
    # m - list of total (population) elements of class i
    # n - list of sampled elements of class i
    #
    # Output:
    # x - multivariate hypergeometric probability for values x_1,x_2,...,x_c
    #
    # Example. From the example given on the web Wikipedia: The Free
    #  Encyclopedia [http://en.wikipedia.org/wiki/Hypergeometric_distribution].
    #  Suppose there are 5 black, 10 white, and 15 red marbles in an urn. You
    #  reach in and randomly select six marbles without replacement. What is
    #  the probability that you pick exactly two of each color?
    #
    #             m = [5,10,15]; n = [2,2,2];
    #
    # Calling on Matlab the function:
    #             x = mhygepdf(m,n)
    #
    # Answer is: (in format long)
    #
    # x =
    #
    #   0.07957559681698
    # Creted by vikash Pandey

    ##### test 'm and must be a list of non-negative and integers.'

    m_int=[isinstance(item, int) for item in m ]
    n_int=[isinstance(item, int) for item in n ]
    if (not all(m_int)):
        print ('M must be a list of non-negative and integers.')
        exit()

    if (not all(n_int)):
        print ('N must be a list of non-negative and integers.')
        exit()
    ### We are going to multivariate compute hypergeometric pdf

    M = np.sum(m)
    N = np.sum(n)
    o = len(m)

    ###
    A = []
    for i in range(o):
        a = np.exp(sc.gammaln(m[i] + 1) - sc.gammaln(n[i] +1) - sc.gammaln(m[i] - n[i] + 1));
        A.append(a)
    ###
    B = np.exp(sc.gammaln(M + 1) - sc.gammaln(N +1) - sc.gammaln(M - N + 1));

    x = max(np.cumprod(A))/B;

    return x
