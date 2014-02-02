DESCRIPTION
-----------

This is the code associated with the following paper:

Strobl EV & Visweswaran S. Markov Blanket Ranking using Kernel-based Conditional Dependence Measures. NIPS Workshop on Causality, 2013.

The algorithm discovers the Markov blanket by backward elimination using a kernel-based conditional dependence measure.

CODE
----

BackElimCD.m - find Markov blanket by backward elimination (more accurate; Algorithm 1 in paper)

ForSelecCD.m - find Markov blanket by forward selection (faster; Algorithm 2)

CreateArtificialMB.m - creates an artificial Markov blanket with linear relationships

demo.m - examples

