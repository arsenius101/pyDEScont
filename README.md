# pyDEScont

## Basics

**pyDEScont** is a Discrete-Event Simulator (DES) for evaluating the performance of buffer-server networks with complex topologies, which are commonly applied to model manufacturing, computer and telecommunications systems with limited buffer space.
In particular, the following basic assumptions are considered: 
- servers are modelled as continuous-time, discrete-state Markov chains. A specific flowrate is assigned to each state of the Markov Chain;
- buffers have finite capacity and are linked to one upstream and one downstream server only; 
- discrete flow in the network (parts to be manufactured, data packets) is considered as continuous.

Flexibility of Markov Chain modelling allows to model servers subject to multiple failure modes, stations made of several servers working in parallel, stages with non-exponential transitions (e.g. Phase-Type), and so on.

## Main features

pyDEScont is an amateur project I started in 2020 in order to learn Python. Its purpose was to convert in Python 3 and improve an old MATLAB script written years ago for my Master of Science Thesis.

These have been the main desired changes from the past:
1) Initialization of system data from an MS Excel file, instead of manual input inside the script itself or an obscure .txt file. **DONE. Excel input/output has been defined in excelIO.py in "lib" folder. A sample input file has been included in the "input_files" folder. Moreover, at startup pyDEScont provides a GUI dialog box to select the system input file.**
2) Parallelization and general script optimization in order to save computation time. **DONE. Although this is far from an accurate benchmark (different machine, different OS), on the same system simulation, on a 4-core processor runtime is about 25-30% of the old (un-parallelized) MATLAB version.**
3) "Recyclable" system initialization module, in order to use the same input for a faster approximated analytical method. **DONE. At the moment, sysdef module, in "lib" folder, has been separated from pyDEScont and will be used for the analytical method (AM) as well.**

## Case studies

Case studies are provided to validate pyDEScont capabilities, comparing its output with results in published papers.
Case studies in sub-folders "1) ..." to "4) ..." are retrieved from **[2]**, while "Split-Merge case" is retrieved from **[1]** (see **"Bibliography"**).
Results are reported in the "SIM_RESULTS" sheet of each Excel file, and can be compared with the relevant tables in the source papers. In the same excels, system and simulation parameters (nr. of runs, duration of each run) are available for transparency.

## Possible improvements

As of December 2020 I was able to create code for an analytical model with the same capabilities of pyDEScont, however the Python algorithm for the two-machine line solver - based on the Tan-Gershwin paper indicated in the **"Bibliography"** section - does not appear to be robust enough when system input data determine very ill-conditioned matrices.
In any case, the code for the two-machine line (aka "building block") is available for in the "lib" folder as "BBsolver.py".
I am willing to share the code for the extended Colledani-Gershwin algorithm I was working on - tentatively named Split-Merge Colledani-Gershwin, or SMCG - if you can help me with the BBsolver. 

## Python version and required libraries

pyDEScont was created using Python 3.8. Required libraries are `numpy`, `scipy`, `pandas`, `openpyxl`.

## License and Contacts

This project is covered by **BSD 3-Clause License**, which allows unlimited redistribution for any purpose as long as its copyright notices and the license's disclaimers of warranty are maintained.
For any further enquiry, my contact e-mail is available in my GitHub profile details.

## Bibliography

**[1]** M. Bolognesi, M. Laudi (2012), _An integrated process-system model for the design and the performance evaluation of recycling systems_, Master of Science Thesis in Industrial Engineering, Politecnico di Milano, Supervisors: T. Tolio, M. Colledani. **Relevant chapter provided in "Case Studies" folder**  

**[2]** M. Colledani, S.B. Gershwin (2013), _A decomposition method for approximate evaluation of continuous flow multi-stage lines with general Markovian machines_, _Annals of Operations Research_ 209(1). [Full-text PDF from MIT website](http://web.mit.edu/manuf-sys/www/oldcell1/papers/colledani-gershwin-anor-2011.pdf)  

**[3]** B. Tan, S.B. Gershwin (2007), _Modeling and Analysis of Markovian Continuous Flow Production Systems with a Finite Buffer: A General Methodology and Applications_, Technical Report ORC-381-07, Massachusetts Institute of Technology Operations Research Center Working Paper Series. [Full-text PDF from MIT website](https://dspace.mit.edu/bitstream/handle/1721.1/37588/ORC-381-07.pdf?sequence=1&isAllowed=y)





