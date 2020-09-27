# pyDEScont

## Basics

**pyDEScont** is a Discrete-Event Simulator (DES) for evaluating the performance of buffer-server networks with complex topologies, which are commonly applied to model manufacturing, computer and telecommunications systems with limited buffer space.
In particular, the following basic assumptions are considered: 
- servers are modelled as continuous-time, discrete-state Markov chains. A specific flowrate is assigned to each state of the Markov Chain;
- buffers have finite capacity and are linked to one upstream and one downstream server only; 
- discrete flow in the network (parts to be manufactured, data packets) is considered as continuous.

Flexibility of Markov Chain modelling allows to model servers subject to multiple failure modes, stations made of several servers working in parallel, stages with non-exponential transitions (e.g. Phase-Type), and so on.

## Main features

pyDEScont is an amateur project started in August 2020. Its purpose was to make a working version in Python 3 starting from a MATLAB script written in 2012 for my Master of Science Thesis in Industrial Engineering.

These have been the main desired changes from the past:
1) Initialization of system data from an MS Excel file, instead of manual input inside the script itself (on an obscure txt file). **DONE. A sample input file has been included in the input_files folder. Moreover, at startup pyDEScont provides a GUI dialog box to select the system input file.**
2) Parallelization and general script optimization in order to save computation time. **DONE. Although this is far from an accurate benchmark (different machine, different OS), on the same system simulation, on a 4-core processor runtime is about 25% of the 2012 version.**
3) "Recyclable" system initialization module, in order to use the same input for a faster approximated analytical method. **Work in progress. At the moment, sysdef module has been separated from pyDEScont and will be used for the analytical method (AM) as well. When pyDEScont is validated, it will be used as a benchmark for the AM accuracy.**

## Next steps

_Modelling capabilities:_
- feasible system layouts include not only lines, but even more complex topologies with multiple branches. However, at the moment it is not capable to model systems with loops (e.g. parts to be reworked due to defects) and that could be implemented in the future;
- it would be interesting to evaluate how well this modelling via Markov Chains can manage non-deterministic flowrates (e.g. exponential or generally-distributed service times).



