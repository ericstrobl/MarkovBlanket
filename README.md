DESCRIPTION
-----------

This is the code associated with the following paper:

Strobl EV & Visweswaran S. Markov Blanket Ranking using Kernel-based Conditional Dependence Measures. NIPS Workshop on Causality, 2013. http://arxiv-web3.library.cornell.edu/abs/1402.0108

The algorithm discovers the Markov blanket by backward elimination or forward selection using a kernel-based conditional dependence measure.

CODE
----

Each file has been designed to run on its own, so you don't need to download everything.

All the functions have a regularization parameters (reg) that we have set to 1E-6 as the default, but we highly recommend that you tune it for each dataset in practice (e.g., try 1E-2,1E-4,1E-6,1E-8).

*Single Markov Blanket*

* BackCD.m - Find Markov blanket by backward elimination. Variables are in "ascending" order. More accurate.

* ForCD.m - Find Markov blanket by forward selection. Variables are in "descending" order. Faster.

*Multiple Markov Blankets*
(top variables are in any reasonable Markov blanket)

* BackCDm.m

* ForCDm.m

*Others*

* CreateArtificialMB.m - Creates an artificial Markov blanket with linear relationships

* demo.m - Examples

