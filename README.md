# GraphPredictiveModels
Use predictive models and ensembles to statistically estimate similarities and differences between cohorts of brain networks. 

Sparse Canonical Correlation Analysis (sCCA) is used as a bi-directional predictive model of brain connectomes (typically structural brain connectomes and functional brain connectomes). The sCCA biconvex criterion implemented in [PMA R toolbox](https://cran.r-project.org/web/packages/PMA/index.html) has been modified based on randomised Lasso principle to allow identification of the most relevant connections based on ensemble/bootstrapping principles. Therefore, relevant structural and functional connections that play an important role to the prediction are identified along with a probability score. 

**Identification.R** projects functional connectivity matrices into an approximate tanget space on the Riemannian manifold, which allows to constrain prediction to Symmetric Positive Definite Matrices (SPD). 


## Related Publications
- F. Deligianni, H. Singh, H.N. Modi, S. Jahani, M. Yucel, A. Darzi, D.R. Leff, G.Z. Yang, 'Expertise and Task Pressure in fNIRS-based brain Connectomes', https://arxiv.org/abs/2001.00114
- F.Deligianni, J. Clayden and G.Z Yang, ['Comparison of Brain Networks Based on Predictive Models of Connectivity'](https://ieeexplore.ieee.org/document/8942010/authors#authors), IEEE BIBE, 2019. (Best Paper Award)
- F. Deligianni, D.W. Carmichael, Gary H. Zhang, C.A. Clark and J.D. Clayden, ['NODDI and tensor-based microstructural indices as predictors of functional connectivity'](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153404), PLoS ONE, 11(4), 2016.
- F. Deligianni, M. Centeno, D.W. Carmichael and J.D. Clayden, ['Relating resting-state fMRI and EEG whole-brain connectomes across frequency bands'](https://www.frontiersin.org/articles/10.3389/fnins.2014.00258/full), Frontiers in Neuroscience, 8(258), 2014. 
