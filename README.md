*Last update: 03/04/14*
DESCRIPTION
-----------

This is code related to kernel-based Markov blanket discovery, either by conditional dependence-based or dependency-based approximation methods.

- *BackCD* and *ForCD* use a conditional dependency-based measure and are from: Strobl EV & Visweswaran S. Markov Blanket Ranking using Kernel-based Conditional Dependence Measures. NIPS Workshop on Causality, 2013. http://arxiv-web3.library.cornell.edu/abs/1402.0108

- *BAHSIC* and *FOHSIC* use a dependency-based measure and are from: Song L, Smola A, Gretton A, Bedo J, and Borgwardt K. Feature Selection via Dependence Maximization. JMLR, 2007.

- *HSIClasso* can be found here: http://www.makotoyamada-ml.com/hsiclasso.html

Tricks of the trade including the Kronecker delta kernel for classification and the copula transform have also been included. Please send me an email if you have any questions/problems: evs17 (at) pitt (dot edu)

CODE
----

Please download the entire package (including the utility functions)

*Main Methods*

- BackCD.m - Find Markov blanket by backward elimination. Variables are in "ascending" order (least to most likely). More accurate.

- ForCD.m - Find Markov blanket by forward selection. Variables are in "descending" order (most to least likely). Faster.

- BAHSIC.m - Backward eliminiation by HSIC as described by Song et al. (2007). Ascending order.

- FOHSIC.m - Forward selection by HSIC as described by Song et al. (2007). Descending order.

*Others*

- CreateArtificialMB.m - Creates an artificial Markov blanket with linear relationships

- demo.m - Examples

