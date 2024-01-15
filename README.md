# Engineering_holography

The code is used to find the optimal graph for two instance of the hyperbolic pentagon code [Pastawski et al.,
JHEP 2015:149 (2015)]:

1) Small instance of 16 qubits (Main_16.py):
   
     * Option 1 (Value='all'): Explore all non-ismorphic graphs using the library gsc* and print the ones with minimal edges.
   
     * Option 2 (Value='ch'): Select the particular set of Hadamard gates to obtain the graph state from the manuscript.  
    
     
2) Bigger instance of 36 qubits (Main_36.py)
   
     * Option 1 (Value='opt'): Found a better graph by applying Hadamard gates in different four pentagons plaquettes a couple of times.
   
     * Option 2 (Value='ch'): Select the particular set of Hadamard gates to obtain the graph state from the manuscript.  

In both cases the program also gives you the files GraphBM.pdf and GraphAM.Pdf. The first shows the layout of the pentagons before the meaurement and the second the resulting graph after the measurement (without optimizing).


*The library gsc is a modified version of https://github.com/sammorley-short/gsc which works for python3.

The preprint of the work is in https://arxiv.org/abs/2209.08954.


