import numpy as np
import Functions as f
import time
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
import LC_explore
from gsc.explore_lc_orbit import apply_qubit_LCs,get_min_edge_reps
import pickle

from gsc.is_lc_equiv import are_lc_equiv

######  Initial data ######

d=200 # edge size

### Building block ###

X=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
Z=np.array([[0,1,0,0,1,1],[1,0,1,0,0,1],[0,1,0,1,0,1],[0,0,1,0,1,1],[1,0,0,1,0,1],[1,1,1,1,1,0]])

### Index contraction ###

Position=[[0,10],[4,12],[6,22],
[9,24],[28,33], [13,40],[36,43],
[31,49],[48,57],[55,64],
[18,16],
[32,1],[44,3],
[63,45],[56,2]] #First you add all the pentagons [i,j] with j>i and then the close loops [i,j] i>j

### Optimizing over Hadamard gates or Optimal Hadamard gates ###

#Value = 'opt'
Value = 'ch'

######  Program ######

### State after measuring ###

GX,GZ,Vph,GraphPosRed=f.StateAM(X,Z,Position, d, figBM=True, layout=True)

### Transforming to a graph state ###

GXwH=GX.copy()
GZwH=GZ.copy()
VphwH=Vph.copy()

GX,GZ,Vph,Halamard=f.GraphTransform(GX,GZ,Vph)

g=f.GraphDraw(GX,GZ)

labels={old_label:new_label for old_label, new_label in enumerate(g.nodes())}
pos0red = {label:new_label for label, new_label in enumerate(GraphPosRed)}

Logq=[0,3,6,10,14,17,21,24,28,31,35] ### logical qubits

arry=np.delete(g.nodes,Logq)

nx.draw_networkx_labels(g, pos0red, labels,font_color='1',font_size=6)
nx.draw_networkx_nodes(g, pos0red, nodelist=Logq, node_color="r",node_size=100)
nx.draw_networkx_nodes(g, pos0red, nodelist=list(arry), node_color="b",node_size=100)
nx.draw_networkx_edges(g,pos0red)

fig=plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig("GraphAM.pdf")
plt.close()

####### Optimization ########

if Value == 'opt':
	#### Optimal Halamard ####

	g_Opt=g.copy()
	
	Opt_H_T=np.zeros(len(g_Opt.nodes()),int)
	
	for rep in range(0,2): ### converge to a local minimum the second time
		
		print("Optimzation %s\n" % rep)
		
		g_Opt, Opt_H= LC_explore.OptHal_36(g_Opt, pos0red, d)
		
		
		nx.draw_networkx_labels(g_Opt, pos0red, labels,font_color='1',font_size=6)
		nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=Logq, node_color="r",node_size=100)
		nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=list(arry), node_color="b",node_size=100)
		nx.draw_networkx_edges(g_Opt,pos0red)

		Graph_name="Graph_%s.pdf" % rep

		fig=plt.gcf()
		fig.set_size_inches(8, 8)
		plt.savefig(Graph_name)

		plt.close()
		
		Opt_H_T=(Opt_H_T+np.array(Opt_H)) % 2 ### summing Hadamard gates
	
	
	### Hadamard gates position  ###

	Opt_H_T=list(Opt_H_T) 

	Hal=np.zeros(len(g.nodes),int)

	for i in range(0,len(g.nodes)):
		if i in Halamard:
			Hal[i]=1

	Opt_H_T=Hal+np.array(Opt_H_T)

	Opt_H_T=list(Opt_H_T)

	Opt_H2_T=[] 

	for k in range(0,len(Opt_H_T)):
		if Opt_H_T[k]==1:
			Opt_H2_T.append(k)

	print("\nHadamard gates=%s" % Opt_H2_T)

	gates=Opt_H2_T.copy()

else:

	#### Particular choose of Halamard ####

	gates=[1, 2, 3, 4, 5, 6, 9, 10, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, 34, 35] # gates from the manuscript
	#gates=[1, 2, 3, 4, 5, 6, 7, 10, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, 32, 35] # gates found in the optimal process
	MX,MZ,Vph=LC_explore.HalamardGates(GXwH,GZwH,VphwH,gates)

	MX,MZ,Vph=f.Triangular(MX,MZ,Vph)
	MX,MZ,Vph=f.CleanMatrix(MX,MZ,Vph)
	g_Opt=f.GraphDraw(MX,MZ)

print(g_Opt)

nx.draw_networkx_labels(g_Opt, pos0red, labels,font_color='1',font_size=6)
nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=Logq, node_color="r",node_size=100)
nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=list(arry), node_color="b",node_size=100)
nx.draw_networkx_edges(g_Opt,pos0red)

fig=plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig("GraphOpt.pdf")
plt.close()

#### Graph with Hadamards in green ###

nx.draw_networkx_labels(g_Opt, pos0red, labels,font_color='1',font_size=6)
nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=Logq, node_color="r",node_size=100)
nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=list(arry), node_color="b",node_size=100)
nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=gates, node_color="g",node_size=100)
nx.draw_networkx_edges(g_Opt,pos0red)

nx.draw_networkx_edges(g_Opt,pos0red)

fig=plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig("GraphOptH.pdf")
plt.close()

### Save information for python3 ###

with open('Graph_36', 'wb') as fp:
	Nodes=list(g_Opt.nodes)
	Edges=list(g_Opt.edges)
	MX,MZ=f.Graphtomatrix(g_Opt)
	pickle.dump([Nodes,Edges,GraphPosRed,labels,MX,MZ], fp)
