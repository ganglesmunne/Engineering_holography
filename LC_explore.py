# Import Python packages
import networkx as nx
# Import local modules
from gsc.explore_lc_orbit import explore_lc_orbit,get_min_edge_reps
import numpy as np
from math import ceil
from gsc.is_lc_equiv import are_lc_equiv
from itertools import permutations, product
import itertools
import Functions as f


### Distance weight ###

def Distance(g,pos0red,d):

	Edg=list(g.edges)

	D=0

	for i in range(0,len(Edg)): 

		cost=np.sqrt((pos0red[Edg[i][0]][0]-pos0red[Edg[i][1]][0])**2+(pos0red[Edg[i][0]][1]-pos0red[Edg[i][1]][1])**2)
		cost=np.exp(ceil(cost/d))
		D=D+cost

	return D

### Permuting nodes in a graph ###

def Permutations(g, Perm, pos0red, d):

	mappingP={old_label:new_label for old_label, new_label in enumerate(Perm)}

	g=nx.relabel_nodes(g, mappingP)

	dist=Distance(g,pos0red,d)

	return g, dist

### Choose one non-ismorphic graph with minimal nodes ###

def Orbit(g,graph_number):

	class_graph = explore_lc_orbit(g)

	g_min=get_min_edge_reps(class_graph)
	
	N=[]
	for obj in g_min:
		N.append(obj)
	
	g=g_min[N[graph_number]]['nx_graph']

	return g

### Exploring all the non-ismorphic graph with minimal nodes ###

def OrbitAll(g):

	# Find the class graph
	class_graph = explore_lc_orbit(g)
	
	g=get_min_edge_reps(class_graph)
	
	#g=int_relabel_graph(g)
	return g


def Optimization_16(g,pos0red,d):


	perm0 = list(permutations([0,1,2,3]))
	perm1 = list(permutations([4,5,6,7]))
	perm2 = list(permutations([8,9,10,11]))
	perm3 = list(permutations([12,13,14,15]))
	
	Nmax=np.inf
	Num=len(perm0)

	for j0,j1,j2,j3 in product(perm0,perm1,perm2,perm3):

		Perm=np.concatenate((j0,j1,j2,j3))
		g_p,dist=Permutations(g, Perm, pos0red , d)
		
		if Nmax>=dist:

			Bin,U=are_lc_equiv(g,g_p)
			
			if Bin == True:
				Nmax=dist
				perm_LC=Perm.copy()
				g_LC=g_p.copy()
				U_LC=U.copy()

	return g_LC,perm_LC,U_LC

	return LCg,LCq,U

def Optimization(g,pos0red,file,d):

	perm0 = list(permutations([0,1,2,3,4,5,6,7,8,9]))
	
	Nmax=np.inf
	Num=len(perm0)

	for j0 in range(0,Num):
		Perm=np.array(perm0[j0])
		g_p,dist=Permutations(g, Perm, pos0red , d)
		
		if Nmax>=dist:

			Bin,U=are_lc_equiv(g,g_p)
			
			if Bin == True:
				Nmax=dist
				perm_LC=Perm.copy()
				g_LC=g_p.copy()
				U_LC=U.copy()
	
	return g_LC,perm_LC,U_LC

	return LCg,LCq,U

def Halam(MX,MZ,Vph,gates):

	Nqub=len(Vph)

	if len(gates)!=0:
		for i in range(0,Nqub):
			if gates[i]==1:
				for j in range(0,Nqub):
					if MX[j][i]==1 and MZ[j][i]==0:
						MX[j][i]=0
						MZ[j][i]=1
					elif MX[j][i]==0 and MZ[j][i]==1:
						MX[j][i]=1
						MZ[j][i]=0
					elif MX[j][i]==1 and MZ[j][i]==1:
						Vph[j]=(Vph[j]+1) %2
	
	return MX,MZ,Vph

def OptHal(g, pos0red, d):

	Nmax=np.inf
	v=np.zeros(len(g.nodes()))
	for gates in itertools.product([0, 1], repeat=len(g.nodes())):

		gatesT=list(gates)
		
		MX,MZ=Halam(g,gatesT)
		MX,MZ,_=f.Triangular(MX,MZ,v)
		MX,MZ,_=f.CleanMatrix(MX,MZ,v)
		Check=f.GraphCheck(MX,MZ)

		if Check==2:

			graph=f.GraphDraw(MX,MZ)

			dist=Distance(graph,pos0red,d)

			if Nmax>=dist:

				Nmax=dist.copy()

				Opt_H=gatesT
				Opt_g=graph

	return Opt_g, Opt_H

def OptHal_36(g,pos0red, d):

	Plaquettes=[[0, 1,2,3, 4,5,6, 7,8,9,10],
	[0, 15,16,17, 29,30,31, 25,26,27,28],
	[0, 4,5,6, 22,23,24, 18,19,20,21],
	[0, 1,2,3, 15,16,17, 11,12,13,14],
	[0, 22,23,24, 29,30,31, 32,33,34,35]]
	
	t=0
	Nmax=np.inf
	Opt_H_F=np.zeros(len(g.nodes),int)
	
	Zeros=np.zeros(len(g.nodes),int)
	
	for Plaq in Plaquettes:
		
		gatesT=list(np.zeros(len(g.nodes),int))
		n_Plaq=len(Plaq)
		
		print("Plaquet %s\n" % t)
		
		for gates in itertools.product([0, 1], repeat=n_Plaq):

			gates=list(gates)

			for k in range(0,n_Plaq):
				gatesT[Plaq[k]]=gates[k]
			
			MX,MZ=f.Graphtomatrix(g)
			MX,MZ,_=Halam(MX,MZ,Zeros,gatesT)
			MX,MZ,_=f.Triangular(MX,MZ,Zeros)
			MX,MZ,_=f.CleanMatrix(MX,MZ,Zeros)
			Check=f.GraphCheck(MX,MZ)

			if Check==2:

				graph=f.GraphDraw(MX,MZ)

				dist=Distance(graph,pos0red,d)

				if Nmax>=dist:
					Nmax=dist
					Opt_H=gatesT.copy()
					Opt_g=graph.copy()
					
		
		Opt_H_F=(Opt_H_F+np.array(Opt_H)) %2
		
		t+=1
		g=Opt_g.copy()

		
	Opt_H_F=list(Opt_H_F)
	return Opt_g, Opt_H_F

def HalamardGates(MX,MZ,Vph,gates):

	Nqub=len(Vph)

	if len(gates)!=0:
		for i in range(0,Nqub):
			if i in gates:
				for j in range(0,Nqub):
					if MX[j][i]==1 and MZ[j][i]==0:
						MX[j][i]=0
						MZ[j][i]=1
					elif MX[j][i]==0 and MZ[j][i]==1:
						MX[j][i]=1
						MZ[j][i]=0
					elif MX[j][i]==1 and MZ[j][i]==1:
						Vph[j]=(Vph[j]+1) %2
	return MX,MZ,Vph


def AppGates(MX,MZ,Vph,gates):

	Nqub=len(Vph)

	if len(gates)!=0:
		for i in range(0,Nqub):
			for k in range(0,len(gates[i])):
				if gates[i][k]=='H':
					
					for j in range(0,Nqub):
						if MX[j][i]==1 and MZ[j][i]==0:
							MX[j][i]=0
							MZ[j][i]=1
						elif MX[j][i]==0 and MZ[j][i]==1:
							MX[j][i]=1
							MZ[j][i]=0
						elif MX[j][i]==1 and MZ[j][i]==1:
							Vph[j]=(Vph[j]+1) %2
		
				elif gates[i][k]=='S':

					for j in range(0,Nqub):
						if MX[j][i]==1 and MZ[j][i]==0:
							MX[j][i]=1
							MZ[j][i]=1
						elif MX[j][i]==1 and MZ[j][i]==1:
							MX[j][i]=1
							MZ[j][i]=0
							Vph[j]=(Vph[j]+1) %2

	return MX,MZ,Vph
