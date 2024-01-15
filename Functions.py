import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


### Check commutation rules ###

def CheckCommutation(Ax,Az,Bx,Bz):
	ChCom=np.matmul(Ax,np.transpose(Bz))-np.matmul(Az,np.transpose(Bx)) % 2
	return ChCom

#### Rotating a buliding block ####

def Rotgraph(d,angle):
	
	#Up#

	posUp=[]

	A=2*np.pi/5
	B=(np.pi-A)/2
	
	x=d*np.cos(A)
	y=d*np.sin(A)
	
	y1=(d/2)*np.tan(B)
	
	r=np.sqrt(y1**2+(d/2)**2)
	
	posUp.append((-d/2, -y1))
	posUp.append((-x-d/2, y-y1))
	posUp.append((0, r))
	posUp.append((x+d/2, y-y1))
	posUp.append((d/2, -y1))
	posUp.append((0, 0))


	#New#
	
	NewPos=[]

	for i in range(0,len(posUp)):
		
		x0=posUp[i][0]
		y0=posUp[i][1]
		
		Central_Point=0
		
		if x0 !=0:

			Beta=np.arctan(abs(y0/x0))
			
			if y0>0 and x0>0:
				Beta=Beta
			
			elif y0>0 and x0<0:
				Beta=np.pi-Beta

			elif y0<0 and x0>0:
				Beta=2*np.pi-Beta
				
			else:
				Beta=np.pi+Beta

		else:
			if y0>0:
				Beta=np.pi/2
				
			elif y0<0:
				Beta=-np.pi/2
			
			else:
				Central_Point=1
		
		if Central_Point==0:
			NewPos.append((r*np.cos(angle+Beta),r*np.sin(angle+Beta)))
		else:
			NewPos.append((0,0))
	return NewPos

### Sum of two generators ###

def Sum2g(g1,g2):

	g=np.zeros((1,len(g1)))

	for i in range(len(g1)):

		g[0][i]=(g1[i]+g2[i]) % 2

	return g

def Phases2g(gx1,gz1,v1,gx2,gz2,v2):

	d=len(gx1)
	
	VT=(v1+v2)% 2
	
	for i in range(0,d):
		
		if gz1[i]*gx2[i]==1:
			VT=(VT+1) % 2
	return VT

### Triangular ###

def Triangular(GX,GZ,Vph):

	Nr=len(GX[0])
	Nc=len(np.transpose(GX)[0])

	for i in range(0,Nr):

		list=[]

		for j in range(0,Nc):

			if j>=i or GX[j][j]==0:

				if GX[j][i]==1 :

					list.append(j)
		if len(list) !=0:

			for k in range(1,len(list)):

				GX[list[k]]=Sum2g(GX[list[k]],GX[list[0]])
				GZ[list[k]]=Sum2g(GZ[list[k]],GZ[list[0]])
				Vph[list[k]]=Phases2g(GX[list[k]],GZ[list[k]],Vph[list[k]],GX[list[0]],GZ[list[0]],Vph[list[0]])
			
			RowX=GX[list[0]]
			RowX2=GX[i]
			
			RowZ=GZ[list[0]]
			RowZ2=GZ[i]
			
			RowV=Vph[list[0]]
			RowV2=Vph[i]

			GX=np.delete(GX,(i),axis=0)
			GX=np.insert(GX,i,RowX,axis=0)

			GX=np.delete(GX,(list[0]),axis=0)
			GX=np.insert(GX,list[0],RowX2,axis=0)

			GZ=np.delete(GZ,(i),axis=0)
			GZ=np.insert(GZ,i,RowZ,axis=0)

			GZ=np.delete(GZ,(list[0]),axis=0)
			GZ=np.insert(GZ,list[0],RowZ2,axis=0)
			
			Vph=np.delete(Vph,(i),axis=0)
			Vph=np.insert(Vph,i,RowV,axis=0)

			Vph=np.delete(Vph,(list[0]),axis=0)
			Vph=np.insert(Vph,list[0],RowV2,axis=0)

		if i==Nc:
			break
	return GX,GZ,Vph

### CleanMatrix ###

def CleanMatrix(GX,GZ,Vph):

	list=[]

	### CleanX ###
	
	for i in range(1,len(np.transpose(GX[0]))):

		if GX[i][i]==1:

			for j in range(1,i+1):

				if GX[i-j][i]==1:

					GX[i-j]=Sum2g(GX[i-j],GX[i])
					GZ[i-j]=Sum2g(GZ[i-j],GZ[i])
					Vph[i-j]=Phases2g(GX[i-j],GZ[i-j],Vph[i-j],GX[i],GZ[i],Vph[i])
	return GX,GZ,Vph

### Graph triangularzation ###

def TriangularZ(GX,GZ,Vph):

	listX=[]

	### Number of 0 in X matrix and Z matrix ###

	for i in range(0,len(np.transpose(GX)[0])):

		if sum(GX[i])==0:
			listX.append(i)

	### Triangularitzation of Z matrix ###

	s=0

	for i in range(0,len(listX)):

		listZ=[]

		for j in range(0,len(listX)):

			if listX[j]>=listX[i] or GZ[listX[j]][listX[j]]==0:
				
				if GZ[listX[j]][listX[i]]==1:

					listZ.append(listX[j])

		if len(listZ) !=0:

			for k in range(1,len(listZ)):

				GX[listZ[k]]=Sum2g(GX[listZ[k]],GX[listZ[0]])
				GZ[listZ[k]]=Sum2g(GZ[listZ[k]],GZ[listZ[0]])
				Vph[listZ[k]]=Phases2g(GX[listZ[k]],GZ[listZ[k]],Vph[listZ[k]],GX[listZ[0]],GZ[listZ[0]],Vph[listZ[0]])
				
			RowX=GX[listZ[0]]
			RowX2=GX[listX[i]]
			
			RowZ=GZ[listZ[0]]
			RowZ2=GZ[listX[i]]
			
			RowV=Vph[listZ[0]]
			RowV2=Vph[listX[i]]

			GX=np.delete(GX,(listX[i]),axis=0)
			GX=np.insert(GX,listX[i],RowX,axis=0)

			GX=np.delete(GX,(listZ[0]),axis=0)
			GX=np.insert(GX,listZ[0],RowX2,axis=0)

			GZ=np.delete(GZ,(listX[i]),axis=0)
			GZ=np.insert(GZ,listX[i],RowZ,axis=0)

			GZ=np.delete(GZ,(listZ[0]),axis=0)
			GZ=np.insert(GZ,listZ[0],RowZ2,axis=0)
			
			Vph=np.delete(Vph,(listX[i]),axis=0)
			Vph=np.insert(Vph,listX[i],RowV,axis=0)

			Vph=np.delete(Vph,(listZ[0]),axis=0)
			Vph=np.insert(Vph,listZ[0],RowV2,axis=0)
			 
	### CleanZ ###

	list=[]

	for i in range(1,len(listX)):

		if GZ[listX[i]][listX[i]]==1:

			for j in range(1,i+1):

				if GZ[listX[i-j]][listX[i]]==1:

					GX[listX[i-j]]=Sum2g(GX[listX[i-j]],GX[listX[i]])
					GZ[listX[i-j]]=Sum2g(GZ[listX[i-j]],GZ[listX[i]])
					Vph[i-j]=Phases2g(GX[listX[i-j]],GZ[listX[i-j]],Vph[listX[i-j]],GX[listX[i]],GZ[listX[i]],Vph[listX[i]])
	return GX,GZ,Vph

### Linear Independent Trace ###

def LI(GX,GZ,Vg):

	GX,GZ,Vg=Triangular(GX,GZ,Vg)
	GX,GZ,Vg=CleanMatrix(GX,GZ,Vg)

	### Trace ###

	list=[]
	
	L=len(np.transpose(GX)[0])

	for i in range(L):
		if sum(GX[i])==0 and sum(GZ[i])==0:
			list.append(i)
	
	if len(list)==0:
		GX,GZ,Vg=TriangularZ2(GX,GZ,Vg)
		list=[]
		for i in range(L):
			
			if sum(GX[i])==0 and sum(GZ[i])==0:
				list.append(i)

	if len(list)==1:

		GX=np.delete(GX,list[0],axis=0)
		GZ=np.delete(GZ,list[0],axis=0)
		Vg=np.delete(Vg,list[0],axis=0)
	else:

		GX=np.delete(GX,(list[0],list[1]),axis=0)
		GZ=np.delete(GZ,(list[0],list[1]),axis=0)
		Vg=np.delete(Vg,(list[0],list[1]),axis=0)
	return GX,GZ,Vg

### Graph merge ###

def GraphMerge(G1,G2,Bell):

	### Check matrix construction ###

	if len(G2) != 0:

		L1=len(G1[0])
		L2=len(G2[0])

		for i in range(2):

			G1[i]=np.transpose(G1[i])
			G1[i]=np.concatenate((G1[i],np.zeros((L2,L1))))
			G1[i]=np.transpose(G1[i])

			G2[i]=np.transpose(G2[i])
			G2[i]=np.concatenate((np.zeros((L1,L2)),G2[i]))
			G2[i]=np.transpose(G2[i])

		GX=np.concatenate((G1[0],G2[0]))
		GZ=np.concatenate((G1[1],G2[1]))
		Vph=np.concatenate((G1[2],G2[2]))

	else:

		GX=G1[0]
		GZ=G1[1]
		Vph=G1[2]

	### Bell measurament ###

	XX=Bell
	ZZ=Bell
	Pos=[]

	for i in range(len(Bell[0])):

		if Bell[0][i]==1:

			Pos.append(i)

	### ZZ ###

	AntiZ=[]


	for i in range(len(ZZ[0])):

		if GX[i][Pos[0]]!=GX[i][Pos[1]]:

			AntiZ.append(i)

	for i in range(len(AntiZ)-1):
 
		GX[AntiZ[i]]=Sum2g(GX[AntiZ[i]],GX[AntiZ[len(AntiZ)-1]])
		GZ[AntiZ[i]]=Sum2g(GZ[AntiZ[i]],GZ[AntiZ[len(AntiZ)-1]])
		Vph[AntiZ[i]]=Phases2g(GX[AntiZ[i]],GZ[AntiZ[i]],Vph[AntiZ[i]],GX[AntiZ[len(AntiZ)-1]],GZ[AntiZ[len(AntiZ)-1]],Vph[AntiZ[len(AntiZ)-1]])

	if len(AntiZ)!=0:

		GX[AntiZ[len(AntiZ)-1]]=np.zeros((1,len(ZZ[0])))
		GZ[AntiZ[len(AntiZ)-1]]=ZZ[0]
		Vph[AntiZ[len(AntiZ)-1]]=0

	### XX ###

	AntiX=[]

	for i in range(len(XX[0])):

		if GZ[i][Pos[0]]!=GZ[i][Pos[1]]:

			AntiX.append(i)    

	for i in range(len(AntiX)-1):

		GX[AntiX[i]]=Sum2g(GX[AntiX[i]],GX[AntiX[len(AntiX)-1]])
		GZ[AntiX[i]]=Sum2g(GZ[AntiX[i]],GZ[AntiX[len(AntiX)-1]])
		Vph[AntiX[i]]=Phases2g(GX[AntiX[i]],GZ[AntiX[i]],Vph[AntiX[i]],GX[AntiX[len(AntiX)-1]],GZ[AntiX[len(AntiX)-1]],Vph[AntiX[len(AntiX)-1]])

	if len(AntiX) != 0:

		GZ[AntiX[len(AntiX)-1]]=np.zeros((1,len(XX[0])))
		GX[AntiX[len(AntiX)-1]]=XX[0]
		Vph[AntiX[len(AntiX)-1]]=0

	### Partial Trace ####

	GX=np.delete(GX,(Pos[0],Pos[1]),axis=1)
	GZ=np.delete(GZ,(Pos[0],Pos[1]),axis=1)
	

	if len(AntiX) != 0 and len(AntiZ) != 0:
		GX=np.delete(GX,([AntiX[len(AntiX)-1]],[AntiZ[len(AntiZ)-1]]),axis=0)
		GZ=np.delete(GZ,([AntiX[len(AntiX)-1]],[AntiZ[len(AntiZ)-1]]),axis=0)
		Vph=np.delete(Vph,([AntiX[len(AntiX)-1]],[AntiZ[len(AntiZ)-1]]),axis=0)
		Com=0

	else:

		#print("Com")
		Com=1

	return GX, GZ, Vph, Com


### Graph Transform ###

def GraphTransform(GX,GZ,Vph):

	Halamard=[]
	list=[]

	for i in range(len(GX[0])):

		if GX[i][i]==0:

			if GZ[i][i]==0:
				
				GX,GZ,Vph=TriangularZ(GX,GZ,Vph)
				break

	for i in range(len(GX[0])):

		if GX[i][i]==0:

			if GZ[i][i]==1:  ### Halamard ###

				Halamard.append(i)

				for j in range(0,len(GX[0])):

					if GX[j][i]==1:

						if GZ[j][i]==0:

							GX[j][i]=0
							GZ[j][i]=1
					
					elif GX[j][i]==1 and GZ[j][i]==1:
						Vph[j]=1

					else:

						if GZ[j][i]==1:

							GX[j][i]=1
							GZ[j][i]=0
					

					
				for j in range(0,len(GX[0])):

					if GX[j][i]==1 and j != i:  

						GX[j]=Sum2g(GX[j],GX[i])
						GZ[j]=Sum2g(GZ[j],GZ[i])
						Vph[j]=Phases2g(GX[j],GZ[j],Vph[j],GX[i],GZ[i],Vph[i])

	return GX,GZ,Vph,Halamard


### Graph without measurment ###

def GraphWmeasurment(G1,G2):
	
	### Check matrix construction ###

	if len(G2) != 0:

		L1=len(G1[0])
		L2=len(G2[0])

		for i in range(2):

			G1[i]=np.transpose(G1[i])
			G1[i]=np.concatenate((G1[i],np.zeros((L2,L1))))
			G1[i]=np.transpose(G1[i])

			G2[i]=np.transpose(G2[i])
			G2[i]=np.concatenate((np.zeros((L1,L2)),G2[i]))
			G2[i]=np.transpose(G2[i])

		GX=np.concatenate((G1[0],G2[0]))
		GZ=np.concatenate((G1[1],G2[1]))

	else:

		GX=G1[0]
		GZ=G1[1]
		
	return GX,GZ


### ###
def TriangularZ2(GX,GZ,Vph):

	### LI of GZ ###

	listX=[]

	### Number of 0 in X matrix ###

	for i in range(0,len(np.transpose(GX)[0])):

		if sum(GX[i])==0:

			listX.append(i)

	### Triangularitzation of Z matrix ###

	Zqub=[]

	s=0
	l=0

	for i in range(len(GZ[0])):

		listZ=[]

		for j in range(0,len(listX)):
			summ=0
			for k in range(j,i):
				summ=summ+GZ[listX[j]][k]

			if j>=i or summ==0:

				if GZ[listX[j]][i]==1:

					listZ.append(listX[j])

		if len(listZ) !=0:

			for j in range(1,len(listZ)):

				GX[listZ[j]]=Sum2g(GX[listZ[j]],GX[listZ[0]])
				GZ[listZ[j]]=Sum2g(GZ[listZ[j]],GZ[listZ[0]])
				Vph[listZ[j]]=Phases2g(GX[listZ[j]],GZ[listZ[j]],Vph[listZ[j]],GX[listZ[0]],GZ[listZ[0]],Vph[listZ[0]])

			RowX=GX[listZ[0]]
			RowX2=GX[listX[s]]
			
			RowZ=GZ[listZ[0]]
			RowZ2=GZ[listX[s]]
			
			RowV=Vph[listZ[0]]
			RowV2=Vph[listX[s]]
			

			GX=np.delete(GX,(listX[s]),axis=0)
			GX=np.insert(GX,listX[s],RowX,axis=0)

			GX=np.delete(GX,(listZ[0]),axis=0)
			GX=np.insert(GX,listZ[0],RowX2,axis=0)

			GZ=np.delete(GZ,(listX[s]),axis=0)
			GZ=np.insert(GZ,listX[s],RowZ,axis=0)

			GZ=np.delete(GZ,(listZ[0]),axis=0)
			GZ=np.insert(GZ,listZ[0],RowZ2,axis=0)
			
			
			Vph=np.delete(Vph,(listX[s]),axis=0)
			Vph=np.insert(Vph,listX[s],RowV,axis=0)

			Vph=np.delete(Vph,(listZ[0]),axis=0)
			Vph=np.insert(Vph,listZ[0],RowV2,axis=0)

			Zqub.append(i)

			s += 1

		if s==len(listX):

			break

	### CleanZ ###

	for i in range(1,len(listX)):

		if GZ[listX[i]][Zqub[i]]==1:
			
			for k in range(1,i+1):
				
				if GZ[listX[i-k]][Zqub[i]]==1:

					GX[listX[i-k]]=Sum2g(GX[listX[i-k]],GX[listX[i]])
					GZ[listX[i-k]]=Sum2g(GZ[listX[i-k]],GZ[listX[i]])
					Vph[listX[i-k]]=Phases2g(GX[listX[i-k]],GZ[listX[i-k]],Vph[listX[i-k]],GX[listX[i]],GZ[listX[i]],Vph[listX[i]])

	return GX,GZ,Zqub,Vph


### Draw a graph from the check matrix ###

def GraphDraw(GX,GZ): 

	g=nx.Graph()

	List_nodes=range(0,len(GZ[0]))

	g.add_nodes_from(List_nodes)
	
	List_edges=[]
	
	for i in range(len(GZ[0])):
		for j in range(len(GZ[0])):

			if GZ[i][j]==1:
				List_edges.append([i,j])

	g.add_edges_from(List_edges)

	return g
 

### Graph to check matrix ###

def Graphtomatrix(g):
	
	n=len(g.nodes())
	
	MX=np.identity(n)
	
	MZ=np.zeros((n,n))
	
	for ed in g.edges():
		MZ[ed[0]][ed[1]]=1
		MZ[ed[1]][ed[0]]=1
	
	return MX,MZ

### Check if it is a graph state ### 

def GraphCheck(GX,GZ):
	Check=0

	### X matrix check ###
	k=0
	for i in range (len(GX[0])):
		for j in range(len(GX[0])):
			
			if i==j and GX[i][j]==1 or i!=j and GX[i][j]==0:
				k+=1

	if k==len(GX[0])**2:
		Check+=1

	k=0

	for i in range (len(GZ[0])):
		ChZ=GZ[i]==np.transpose(GZ)[i]
		
		if ChZ.all()==True :
			k+=1

	if k==len(GX[0]):
		Check+=1

	return Check


### Graph before measuring ###

def GraphBM(X,Z,Position,d,fig,layout):

	Xg0=X
	Zg0=Z

	####### Qubit Position ########

	QuP=range(0,2*len(X))
	TQuP=Position[0]

	for i in range(1,len(Position)):

		if Position[i][1]>Position[i][0]:

			QuP=np.append(QuP,range(len(QuP),len(QuP)+len(X)))

		TQuP=np.append(TQuP,Position[i])

	####### Graphic without measurment ########

	for j in range(len(Position)):
		
		if Position[j][1]>Position[j][0]:

			Gr0=[X,Z]

		else:

			Gr0=[]

		G0=[Xg0,Zg0]
		Xg0,Zg0=GraphWmeasurment(G0,Gr0)


	####### Graph position ########
	
	g0=GraphDraw(Xg0,Zg0)
	
	if layout == True:
		
		Rad=8*d/10

		A=2*np.pi/5
		B=(np.pi-A)/2

		h=(d/2)*np.tan(B)
		y=d*np.sin(A)
		x=d*np.cos(A)
		C=np.pi/2+np.arctan(abs((y-h)/(x+d/2)))
		l=(d/2)/np.cos(B)

		SpaceP=[(0,0,0,d),(2*B,-3*d/2,-3*h,d),(-2*B,3*d/2,-3*h,d),(np.pi, 0 ,-6*h, d),
		(3*A/2,-6*h*np.sin(A),-6*h*np.cos(A),d),(-C,-3*(x+d/2),3*(y-h),d), (-3*A/2,6*h*np.sin(A),-6*h*np.cos(A),d),(C,3*(x+d/2),3*(y-h),d),
		(3*A/2,-6*h*np.sin(A/2),6*h*np.cos(A/2),d),(np.pi, 0 , 3*l, d),(-3*A/2,6*h*np.sin(A/2),6*h*np.cos(A/2),d)]
		 
		
	else:
		
		Rad=7*d/8
		SpaceP=[(0,0,0,d),(-np.pi/2,2*Rad,-2*Rad,d),(np.pi, 0 , -4*Rad, d),(np.pi/2,-2*Rad,-2*Rad,d)]


	GraphPos=[]

	for i in range(len(SpaceP)):

		RotPos=Rotgraph(SpaceP[i][3],SpaceP[i][0])

		for j in range(len(X)):

			GraphPos.append((RotPos[j][0]+SpaceP[i][1],RotPos[j][1]+SpaceP[i][2]))


	if fig == True:
		
		pos0 = {label:new_label for label, new_label in enumerate(GraphPos)}
		
		Pos_array=np.array(Position)

		Pos_array=list(Pos_array.ravel())
		Pos_array.sort()

		arry=np.delete(g0.nodes,Pos_array)
		labels={old_label:new_label for old_label, new_label in enumerate(g0.nodes())}

		nx.draw_networkx_labels(g0, pos0, labels,font_color='1',font_size=6)
		nx.draw_networkx_nodes(g0, pos0, nodelist=Pos_array, node_color="g",node_size=100)
		nx.draw_networkx_nodes(g0, pos0, nodelist=list(arry), node_color="b",node_size=100)
		nx.draw_networkx_edges(g0,pos0)

		gcon=nx.Graph()

		gcon.add_nodes_from(g0.nodes())

		List_edges=Position

		gcon.add_edges_from(List_edges)

		nx.draw_networkx_edges(gcon,pos0,edge_color="r")

		fig=plt.gcf()
		fig.set_size_inches(8, 8)
		plt.savefig("GraphBM.pdf")
		plt.close()

	return GraphPos,TQuP


### Stabilizer state after measuring ###

def StateAM(X,Z,Position,d,figBM,layout): 
	
	
	#FigBM add as a pdf the figure of the graph before measuring
	#Layout True gives you the hyperbolic layout for 36 qubits. False give you the layout for 16 qubits.

	
	GraphPos,TQuP=GraphBM(X,Z,Position,d,figBM, layout) 


	####### Merge ########

	Vr=np.zeros(len(X))

	Xg=X
	Zg=Z
	Vg=Vr

	for i in range(len(Position)):

		PTr=[]

		if Position[i][1]>Position[i][0]:

			Gr=[X,Z,Vr]

			n=2*len(X[0])+i*len(X[0])

		else:

			Gr=[]

		Bell=np.zeros((1,n))

		for j in range(i+1):

			Bell[0][Position[j][0]]=1
			Bell[0][Position[j][1]]=1

		for k in range(i):

			PTr.append(Position[k][0]) 
			PTr.append(Position[k][1])
		
		Bell=np.delete(Bell,PTr,axis=1)
		
		G=[Xg,Zg,Vg]

		Xg,Zg,Vg,Com=GraphMerge(G,Gr,Bell)

		if Com==1:
			Xg,Zg,Vg=LI(Xg,Zg,Vg)
	
	### Clean generators ###
	
	Xg,Zg,Vg=Triangular(Xg,Zg,Vg)

	Xg,Zg,Vg=CleanMatrix(Xg,Zg,Vg)
	
	### Graph position after the contraction ###

	GraphPosRed=GraphPos[:]

	TQuP=list(TQuP)
	TQuP.sort()

	for i in range(len(TQuP)):

		j=(len(TQuP)-1)-i
		GraphPosRed.pop(TQuP[j])
	
	return Xg,Zg,Vg,GraphPosRed
