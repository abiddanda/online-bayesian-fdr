## February 15, 2016

Back to the context of correlated tests in GWAS I came across a pretty interesting paper that maybe we could borrow from / explore (possibly even for the online testing stuff as well). 

http://biorxiv.org/content/biorxiv/early/2016/02/09/039099.full.pdf

Essentially they set up a hmm with hidden states being a true association or not and transition probabilities based off LD (more specifically related to probability of recombination events within a given distance between markers). The emission probabilities were based off of whatever test-statistic used for the association in both the case-control discrete context and quantitative trait mapping. What is kind of interesting is that one can summarize this info with posterior probabilities of a marker being in a particular state rather than p-values and thus they argue avoids some of the multiple testing problems. 

I’m thinking we can kind of adapt this to the bayesian mixture model setup for fdr 

* transition probabilities based off some known correlation for whatever problem at hand (potentially time dependent correlation for online testing / LD for GWAS)
* emission probabilities based off the normal component for signals and the normal component for nulls. We can think about this in the context of beta_hats estimates from gwas and have the signal component be a mixture itself reflecting positive and negative modes. Also use std error summary statistics in variance normal components
* We can also think about how exactly to estimate the proportion of nulls in this context (could be reflected in the emission probabilities)

This would allow both the mixture model to be estimated as well as give a sort of measure of testing via the posteriors (i.e. forward-backwards in the hmm)

Estimation > testing! 

