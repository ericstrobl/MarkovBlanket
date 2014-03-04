*Last update: 03/04/14*
DESCRIPTION
-----------

This is the code associated with the following paper:

Strobl EV & Visweswaran S. Markov Blanket Ranking using Kernel-based Conditional Dependence Measures. NIPS Workshop on Causality, 2013. http://arxiv-web3.library.cornell.edu/abs/1402.0108

The algorithms discover the Markov blanket by backward elimination or forward selection using a kernel-based conditional dependence measure.

CODE
----

Please download the entire package (including the utility functions)

*Main Methods*

* BackCD.m - Find Markov blanket by backward elimination. Variables are in "ascending" order (least to most likely). More accurate.

* ForCD.m - Find Markov blanket by forward selection. Variables are in "descending" order (most to least likely). Faster.

* BAHSIC.m - Backward eliminiation by HSIC as described by Song et al. (2007). Ascending order.

* ForCD.m - Forward selection by HSIC as described by Song et al. (2007). Descending order.

*Others*

* CreateArtificialMB.m - Creates an artificial Markov blanket with linear relationships

* demo.m - Examples

