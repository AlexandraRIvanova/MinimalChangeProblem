
# Minimal change parameter synthesis problem

There are **four** arguments, that are obligatory:
- path to a model file - must be a pMC
- path to property file or property in string fortmat
- initial instantiation in the form as a region, e.g. 0.519\<=p<=0.537,-0.013<=q<=0.005
- epsilon - double

We provide the following options. A user can secify the values of:
- depth limit **x** - in how many regions a region will be divided. (**x** from 1 to *none* - int)
- refinement factor **y** - how many percent of the initial region can be left *unknown*. (**y** between 0 and 1 - double)
- number of random **z** - how many instantiations will be extracted in the end to compare with eachother to find the one with minimal distance (**z** from 1 -int)
 
If not specified the default values are:
- depth limit - *none*
- refinement factor - 0.05
- number of random - 3

We added an option regarding pMCs of BNs: \
- faster **true** - after the second iteration of PLA, tries the whole parameter space to deem the constraint infeasible faster. That options influences the computation time.
