# SGD
This repository contains simulations and empirical experiments of SGD.

'app.R' is an R Shiny app that visualizes the descent of a few different kinds of SGD algorithms (Basic fixed-step, basic diminishing-step, momentum, and ADAM). I created this for a presentation in the Stochastic Optimization seminar at Purdue.

'LogistRegSGD.R' allows for running (and comparing metrics for) different SGD algorithms for L2 regularized logistic regression. Adaptive algorithms of my current research, some standard SGD algorithms, and the Newton-CG algorithm (Nocedal et al. 2018) are included.
