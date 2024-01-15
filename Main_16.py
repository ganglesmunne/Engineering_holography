import numpy as np
import Functions as f
import networkx as nx
import matplotlib.pyplot as plt
import LC_explore
import pickle

######  Initial data ######

d=200 # edge size

### Building block ###

X=np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
Z=np.array([[0,1,0,0,1,1],[1,0,1,0,0,1],[0,1,0,1,0,1],[0,0,1,0,1,1],[1,0,0,1,0,1],[1,1,1,1,1,0]])

Position=[[4,6],[10,12],[16,18],[22,0]] #First you add all the pentagons [i,j] with j>i and then the close loops [i,j] i>j

### Explore all non-ismorphic orbits with minimal edges or Optimal graph from a specific orbit (permutations) or  Optimal Hadamard gates ###

Value = 'all'
#Value = 'ch'

### State after measuring ###

GX,GZ,Vph,GraphPosRed=f.StateAM(X,Z,Position,d,figBM=True, layout=False)

### Transforming to a graph state ###

GXwH=GX.copy()
GZwH=GZ.copy()
VphwH=Vph.copy()

GX,GZ,Vph,Halamard=f.GraphTransform(GX,GZ,Vph)

g=f.GraphDraw(GX,GZ)

labels={old_label:new_label for old_label, new_label in enumerate(g.nodes())}
pos0red = {label:new_label for label, new_label in enumerate(GraphPosRed)}

Logq=[3,7,11,15] ### logical qubits

arry=np.delete(g.nodes,Logq)

nx.draw_networkx_labels(g, pos0red, labels,font_color='1',font_size=8)
nx.draw_networkx_nodes(g, pos0red, nodelist=Logq, node_color="r",node_size=150)
nx.draw_networkx_nodes(g, pos0red, nodelist=list(arry), node_color="b",node_size=150)
nx.draw_networkx_edges(g,pos0red)

fig=plt.gcf()
fig.set_size_inches(8, 8)
plt.savefig("GraphAM.pdf")
plt.close()


####### Optimization ########

if Value == 'all':

	##### Representation all graph non-isomorphic with minimal number of edges ######
	
	graphOrbit=LC_explore.OrbitAll(g)

	List_nodes_Total=range(0,len(g.nodes()))

	for obj in graphOrbit: 

		g=graphOrbit[obj]['nx_graph']

		nx.draw_networkx_labels(g, pos0red, labels,font_color='1',font_size=8)
		nx.draw_networkx_nodes(g, pos0red, nodelist=Logq, node_color="r",node_size=150)
		nx.draw_networkx_nodes(g, pos0red, nodelist=list(arry), node_color="b",node_size=150)
		nx.draw_networkx_edges(g,pos0red)
		print("n%s: %s\n" % (obj,g))
		fig=plt.gcf()
		fig.set_size_inches(8, 8)
		plt.savefig("Graph_%s.pdf" % obj)
		plt.close()

else:

	##### Optimal Graph ####

	gates=[2,6,10,14] ### gates from the manuscript

	MpX,MpZ,Vph=LC_explore.HalamardGates(GXwH,GZwH,Vph,gates)
	MpX,MpZ,Vph=f.Triangular(MpX,MpZ,Vph)
	MpX,MpZ,Vph=f.CleanMatrix(MpX,MpZ,Vph)

	g_Opt=f.GraphDraw(MpX,MpZ)

	print(g_Opt)

	Nodes=list(g_Opt.nodes)
	Edges=list(g_Opt.edges)

	labels={old_label:new_label for old_label, new_label in enumerate(g.nodes())}
	pos0red = {label:new_label for label, new_label in enumerate(GraphPosRed)}

	arry=np.delete(g_Opt.nodes,Logq)

	nx.draw_networkx_labels(g_Opt, pos0red, labels,font_color='1',font_size=8)
	nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=Logq, node_color="r",node_size=150)
	nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=list(arry), node_color="b",node_size=150)
	nx.draw_networkx_edges(g_Opt,pos0red)

	fig=plt.gcf()
	fig.set_size_inches(8, 8)
	plt.savefig("GraphOpt.pdf")
	plt.close()
	

	#### Graph with Hadamards in green ###

	nx.draw_networkx_labels(g_Opt, pos0red, labels,font_color='1',font_size=8)
	nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=Logq, node_color="r",node_size=150)
	nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=list(arry), node_color="b",node_size=150)
	nx.draw_networkx_nodes(g_Opt, pos0red, nodelist=gates, node_color="g",node_size=150)
	nx.draw_networkx_edges(g_Opt,pos0red)

	nx.draw_networkx_edges(g_Opt,pos0red)

	fig=plt.gcf()
	fig.set_size_inches(8, 8)
	plt.savefig("GraphOptH.pdf")
	plt.close()



