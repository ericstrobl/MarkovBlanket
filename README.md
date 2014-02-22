*Last update: 02/07/14*
DESCRIPTION
-----------

This is the code associated with the following paper:

Strobl EV & Visweswaran S. Markov Blanket Ranking using Kernel-based Conditional Dependence Measures. NIPS Workshop on Causality, 2013. http://arxiv-web3.library.cornell.edu/abs/1402.0108

The algorithms discover the Markov blanket by backward elimination or forward selection using a kernel-based conditional dependence measure.

CODE
----

Each file has been designed to run on its own, so you don't need to download everything.

*Main Methods*

* BackCD.m - Find Markov blanket by backward elimination. Variables are in "ascending" order. More accurate.

* ForCD.m - Find Markov blanket by forward selection. Variables are in "descending" order. Faster.


*Others*

* CreateArtificialMB.m - Creates an artificial Markov blanket with linear relationships

* demo.m - Examples

