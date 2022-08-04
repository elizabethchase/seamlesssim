This is the repository for the seamlesssim package, which accompanies Boonstra, Braun, and Chase's paper on seamless trial design ([Clinical Trials 2021](https://journals.sagepub.com/doi/full/10.1177/1740774520981939)). Interested users should examine that paper, along with the manual and vignette posted here, for examples and explanations for how to use the package. 

Users should first check to make sure they have Stan and the rstan package installed. Installing these packages can be challenging. 

If you DON'T use Mac Catalina: We recommend that you go to https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and follow the instructions given there. 

If you DO use Mac Catalina: After Mac updated to Catalina, installing Stan has become something of a chore. We'd recommend using the installation instructions given here: https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-a-Mac If you're still having trouble, you might try the following help pages:

- https://github.com/jeroen/V8/issues/59 

- https://discourse.mc-stan.org/t/dealing-with-catalina/11285

- https://discourse.mc-stan.org/t/dealing-with-catalina-ii/11802

- https://discourse.mc-stan.org/t/dealing-with-catalina-iii/12731

- https://discourse.mc-stan.org/t/dealing-with-catalina-iv/13502


After rstan has been successfully installed, our package can be installed by running the command `install.packages("devtools")` followed by `devtools::install_github("elizabethchase/seamlesssim")`. It will print many error messages during installation--these are not actually errors and not cause for concern. If at the end of installation you see `DONE (seamlesssim)` (as opposed to `non-zero exit status` or `Failed to install 'seamlesssim' from GitHub`), you should be all set. 

If you are still having trouble, we recommend checking that you have updated to the latest versions of both R and RStudio and then attempt reinstalling seamlesssim. 
