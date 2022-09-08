# minimal change problem
This is the implementation of the approach, which is reported in my Bachler's thesis.
The code is provided with the purpose of future experiments on possible improvements.
It is built upon STORM - a model checker provided by the chair MOVES of the RWTH university.

Running requires the following arguments:
- path to the model(DTMC)
- property for the model repair
- epsilon - a bound for the distance between the repaired and original model
- initial instantiation in form of region "0<=p<=0.9,0.1<=q<=0.7"

There are several options for running:
 --numberOfRandom (int) - the number of instantiations, that the method will compare for the minimal one
 --refThreshold (double < 1) - the refinement threshold given to PLA in STORM, indicates how much of the region can be left Unknown
 --depthLimit - the depth limit given to PLA in STORM, indicates the maximum number of times a region will be divided.
