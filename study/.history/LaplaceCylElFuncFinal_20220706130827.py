""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""
   
# LaplaceLine.py:  Solve Laplace's eqtn, 3D matplot, close shell to quit

from mpi4py import MPI
import matplotlib.pylab as p;
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.sparse import csc_matrix
import time
import sys
import gc
import pandas as pd

def plotPot(Val, numx, numy, dt):                                         # V(x, y) 
	x = range(0, numx, dt);  y = range(0, numy, dt)
	X, Y = p.meshgrid(x,y)                 
	Z = Val[X,Y]                        
	fig = p.figure()                                      # Create figure
	ax = Axes3D(fig)                                      # Plot axes
	ax.plot_wireframe(X, Y, Z, color = 'r')               # Red wireframe
	ax.set_xlabel('X')                                     
	ax.set_ylabel('Y')
	ax.set_zlabel('Potential')
	p.show()                                              # Show fig

def InitElec(Filename, Den, Nx, RFV):

	df = pd.read_csv(Filename, index_col=0)
#	print("df = ", df)
	for index, row in df.iterrows():
#		print("index, row = ", index, row)
		Den.extend([[index, float(row['Vol'])+RFV*float(row['RF']), float(row['X1'])/Nx, float(row['Y1'])/Nx, float(row['X2'])/Nx, float(row['Y2'])/Nx]])
#	print("Den = ", Den)
	del df
	gc.collect()

def Electrode(Pot, Elec, Den, numx, numy, offsetx, deltax, t):

#	print(Den)
	for e in Den:
		i, V, x1, y1, x2, y2  = e
#		if myrank == 0:
#			print("e, x1, y1 ... = ", e, x1, y1, x2, y2, i, V)
			
		x1 = int(x1*numx)	#Always lower left
		y1 = int(y1*numx)
		x2 = int(x2*numx)	#Always upper right
		y2 = int(y2*numx)
#		print("x1, x2, y1, y2 = ", x1, x2, y1, y2)
		dx = 0
		if x2-x1 == 0:
			dx = 1
		dy = 0
		if y2-y1 == 0:
			dy = 1

		if x1 < offsetx:
			if x2 < offsetx:
				continue
			elif offsetx <= x2 < offsetx + deltax:
				x1 = 0
				x2 = x2 - offsetx + dx
			else:
				x1 = 0
				x2 = deltax
		elif offsetx <= x1 < offsetx + deltax:
			if offsetx <= x2 < offsetx + deltax:
				x1 -= offsetx
				x2 = x2 -offsetx + dx
			else:
				x1 -= offsetx
				x2 = deltax
		else:
			continue

		if y1 < 0:
			if y2 < 0:
				continue
			elif 0 <= y2 < numy:
				y1 = 0
				y2 = y2 + dy
			else:
				y1 = 0
				y2 = numy
		elif 0 <= y1 < numy:
			if 0 <= y2 < numy:
				y2 = y2 + dy
			else:
				y2 = numy
		else:
			continue

		Elec[x1:x2, y1:y2] = i
		Pot[x1:x2, y1:y2] = V

#		if myrank == 0:
#			print("Elec = ", Elec)
#			print("Pot = ", Pot)

def SplitArray(Nsize, size, Gsx, Olap):

	Nsize.clear()
	if size == 1:
		Nsize.append(Gsx)

	elif 1 < size < 5:
		half_Nlmax = ((Gsx+2*Olap+(size-2)*2*Olap) /(2* size))

		half_delta = half_Nlmax - int(half_Nlmax)
		if half_delta < 0.1:
			Nsize.append((int(half_Nlmax))*2)
			sum = (int(half_Nlmax))*2
			delta = 0.5
		else:
			Nsize.append((int(half_Nlmax)+1)*2)
			sum = (int(half_Nlmax)+1)*2
			delta = 0.5 + half_delta
			if half_delta > 0:
				while(delta<1):
					delta += half_delta
			delta -= 1
		for i in range(size-2):
			delta += half_delta
	#		print("delta, sum ", delta, sum)
			if( delta > 1):
				delta -= 1
				Nsize.append((int(half_Nlmax)+1)*2)
				sum += (int(half_Nlmax)+1)*2
			else:
				Nsize.append((int(half_Nlmax))*2)
				sum += int(half_Nlmax)*2
		
		Nsize.append(int((Gsx+2*Olap+(size-2)*2*Olap) - sum))

	elif size >= 5:
		half_Nlmax = ((Gsx+2*Olap+(size-2)*2*Olap) /(2* size))
		Nsize.append((int(half_Nlmax)+1)*2)
		Nsize.append((int(half_Nlmax)+1)*2)
		sum = 0
		half_Nlmax = (Gsx+2*Olap+(size-2)*2*Olap - 4*int(half_Nlmax+1)) /(2* (size-2))
		half_delta = half_Nlmax - int(half_Nlmax)

		delta = 0.5
		for i in range(size-2):
			delta += half_delta
	#		print("delta, sum ", delta, sum)
			if( delta > 1):
				delta -= 1
				Nsize.insert(1, (int(half_Nlmax)+1)*2)
				sum += (int(half_Nlmax)+1)*2
			else:
				Nsize.insert(1, (int(half_Nlmax))*2)
				sum += int(half_Nlmax)*2

	print("Nsize = ", Nsize)

def SetCSC(gsx, gsy, type):
	"""
	type 0:size=1, type 1:myrank=0(first), type 2:myrank=size-1(last), type 3:between
	Create sparse matrix for SOR method
	"""
#	RA.clear()
#	BA.clear()

	Rrow_index = []
	Rcol_index = []
	Rvalue = []
	Brow_index = []
	Bcol_index = []
	Bvalue = []

	for i in range(gsx):
		for j in range(gsy):
			k = i*gsy+j
			if Elec[i,j] != 0:
				Rrow_index.append(k)
				Rcol_index.append(k)
				Rvalue.append(1.0)
				Brow_index.append(k)
				Bcol_index.append(k)
				Bvalue.append(1.0)
			elif (i + j) % 2 == 0:
				Brow_index.append(k)
				Bcol_index.append(k)
				Bvalue.append(1.0)
				if i == 0: #Left
					if type == 0 or type == 1: #Actual Left
						if j == 0:
							Rrow_index.append(k)
							Rcol_index.append((i + 1)*gsy+j)
							Rvalue.append(0.2)

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j + 1)
							Rvalue.append(0.8)

						elif j == gsy-1:
							Rrow_index.append(k)
							Rcol_index.append((i + 1)*gsy+j)
							Rvalue.append(2*j/(4*j-1))

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j - 1)
							Rvalue.append((2*j-1)/(4*j-1))
						else:
							Rrow_index.append(k)
							Rcol_index.append((i + 1)*gsy+j)
							Rvalue.append(0.33333333)

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j - 1)
							Rvalue.append(0.16666666*(2-1/j))

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j + 1)
							Rvalue.append(0.16666666*(2+1/j))
					else:
						Rrow_index.append(k)
						Rcol_index.append(k)
						Rvalue.append(1.0)
				elif i == gsx-1: #Right
					if type == 0 or type == 2: #Actual Right
						if j == 0:
							Rrow_index.append(k)
							Rcol_index.append((i - 1)*gsy+j)
							Rvalue.append(0.2)

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j + 1)
							Rvalue.append(0.8)
						elif j == gsy-1:
							Rrow_index.append(k)
							Rcol_index.append((i - 1)*gsy+j)
							Rvalue.append(2*j/(4*j-1))

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j - 1)
							Rvalue.append((2*j-1)/(4*j-1))
						else:
							Rrow_index.append(k)
							Rcol_index.append((i - 1)*gsy+j)
							Rvalue.append(0.33333333)

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j - 1)
							Rvalue.append(0.16666666*(2-1/j))

							Rrow_index.append(k)
							Rcol_index.append(i*gsy+j + 1)
							Rvalue.append(0.16666666*(2+1/j))
					else:
						Rrow_index.append(k)
						Rcol_index.append(k)
						Rvalue.append(1.0)

				elif j == 0: #on Axis
					Rrow_index.append(k)
					Rcol_index.append((i - 1)*gsy+j)
					Rvalue.append(0.1666666)

					Rrow_index.append(k)
					Rcol_index.append((i + 1)*gsy+j)
					Rvalue.append(0.1666666)

					Rrow_index.append(k)
					Rcol_index.append(i*gsy+j + 1)
					Rvalue.append(0.6666666)
				elif j == gsy-1:
					Rrow_index.append(k)
					Rcol_index.append((i - 1)*gsy+j)
					Rvalue.append(2*j/(6*j-1))

					Rrow_index.append(k)
					Rcol_index.append((i + 1)*gsy+j)
					Rvalue.append(2*j/(6*j-1))

					Rrow_index.append(k)
					Rcol_index.append(i*gsy+j - 1)
					Rvalue.append((2*j-1)/(6*j-1))
				else:
					Rrow_index.append(k)
					Rcol_index.append((i - 1)*gsy+j)
					Rvalue.append(0.25)

					Rrow_index.append(k)
					Rcol_index.append((i + 1)*gsy+j)
					Rvalue.append(0.25)

					Rrow_index.append(k)
					Rcol_index.append(i*gsy+j - 1)
					Rvalue.append(0.125*(2-1/j))

					Rrow_index.append(k)
					Rcol_index.append(i*gsy+j + 1)
					Rvalue.append(0.125*(2+1/j))
			else:
				Rrow_index.append(k)
				Rcol_index.append(k)
				Rvalue.append(1.0)
				if i == 0: #Left
					if type == 0 or type == 1: #Actual Left
						if j == 0:
							Brow_index.append(k)
							Bcol_index.append((i + 1)*gsy+j)
							Bvalue.append(0.2)

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j + 1)
							Bvalue.append(0.8)

						elif j == gsy-1:
							Brow_index.append(k)
							Bcol_index.append((i + 1)*gsy+j)
							Bvalue.append(2*j/(4*j-1))

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j - 1)
							Bvalue.append((2*j-1)/(4*j-1))
						else:
							Brow_index.append(k)
							Bcol_index.append((i + 1)*gsy+j)
							Bvalue.append(0.33333333)

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j - 1)
							Bvalue.append(0.16666666*(2-1/j))

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j + 1)
							Bvalue.append(0.16666666*(2+1/j))
					else:
						Brow_index.append(k)
						Bcol_index.append(k)
						Bvalue.append(1.0)
				elif i == gsx-1: #Right
					if type == 0 or type == 2: #Actual Right
						if j == 0:
							Brow_index.append(k)
							Bcol_index.append((i - 1)*gsy+j)
							Bvalue.append(0.2)

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j + 1)
							Bvalue.append(0.8)
						elif j == gsy-1:
							Brow_index.append(k)
							Bcol_index.append((i - 1)*gsy+j)
							Bvalue.append(2*j/(4*j-1))

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j - 1)
							Bvalue.append((2*j-1)/(4*j-1))
						else:
							Brow_index.append(k)
							Bcol_index.append((i - 1)*gsy+j)
							Bvalue.append(0.33333333)

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j - 1)
							Bvalue.append(0.16666666*(2-1/j))

							Brow_index.append(k)
							Bcol_index.append(i*gsy+j + 1)
							Bvalue.append(0.16666666*(2+1/j))
					else:
						Brow_index.append(k)
						Bcol_index.append(k)
						Bvalue.append(1.0)
				elif j == 0: #on Axis
					Brow_index.append(k)
					Bcol_index.append((i - 1)*gsy+j)
					Bvalue.append(0.1666666)

					Brow_index.append(k)
					Bcol_index.append((i + 1)*gsy+j)
					Bvalue.append(0.1666666)

					Brow_index.append(k)
					Bcol_index.append(i*gsy+j + 1)
					Bvalue.append(0.6666666)
				elif j == gsy-1:
					Brow_index.append(k)
					Bcol_index.append((i - 1)*gsy+j)
					Bvalue.append(2*j/(6*j-1))

					Brow_index.append(k)
					Bcol_index.append((i + 1)*gsy+j)
					Bvalue.append(2*j/(6*j-1))

					Brow_index.append(k)
					Bcol_index.append(i*gsy+j - 1)
					Bvalue.append((2*j-1)/(6*j-1))
				else:
					Brow_index.append(k)
					Bcol_index.append((i - 1)*gsy+j)
					Bvalue.append(0.25)

					Brow_index.append(k)
					Bcol_index.append((i + 1)*gsy+j)
					Bvalue.append(0.25)

					Brow_index.append(k)
					Bcol_index.append(i*gsy+j - 1)
					Bvalue.append(0.125*(2-1/j))

					Brow_index.append(k)
					Bcol_index.append(i*gsy+j + 1)
					Bvalue.append(0.125*(2+1/j))


	Rs = (gsx*gsy, gsx*gsy)
	Bs = (gsx*gsy, gsx*gsy)

	RA = csc_matrix((Rvalue, (Rrow_index, Rcol_index)), Rs, dtype=np.float32)
	BA = csc_matrix((Bvalue, (Brow_index, Bcol_index)), Bs, dtype=np.float32)

	return(RA, BA)

def Iteration(V, RA, BA, Nsize, comm, size, myrank, gsy, niter, diter, olap):
	buf = np.empty((gsy*olap, 1), np.float32)
	delta = 5
	gsx = Nsize[myrank]
	for iter in range(0, niter, diter):
		Vin = V.copy()
		if myrank == 0:
			if myrank == 0 and iter % 5 < diter:
				print(iter)
				print("Delta, acc = ", delta, acc)

		# The actual iteration
		for initer in range(diter):
			V = RA.dot(V)

			Vd = V - Vin
			Rdelta = np.max(abs(Vd))

			V = Vin + Vd*(1+acc)

			Vin = V.copy()
			V = BA.dot(V)
			Vd = V - Vin
			Bdelta = np.max(abs(Vd))

			V = Vin + Vd*(1+acc)
			Vin = V.copy()


	#		acc = Maxacc

		'''
		if iter < 1:
			acc = 0.4
		elif iter < 10:
			acc = 0.67
		else:
			Hist = 0.7; Fmax = 0.9
			acc = acc*Hist + (1-Hist)*(Fmax - (Fmax-0.4)*(1/(1+delta/(2*0.005))))
		'''

		if Rdelta > Bdelta:
			delta = Rdelta
		else:
			delta = Bdelta
		teldelta[0] = delta
		comm.Reduce(teldelta, alldelta, op=MPI.MAX, root=0)
		delta = alldelta[0]

		if size != 1:
			if myrank == 0:
				reqsend0 = comm.Isend(V[(gsx-2*olap)*gsy:(gsx-olap)*gsy], dest=1, tag=0)
				reqrecv0 = comm.Irecv(buf, source=1, tag=1)
				reqsend0.Wait()
				reqrecv0.Wait()
				V[(gsx-olap)*gsy:(gsx)*gsy] = buf
			elif myrank == size-1:
				reqrecv0 = comm.Irecv(buf, source=size-2, tag=size-2)
				reqsend0 = comm.Isend(V[olap*gsy:(2*olap)*gsy], dest=size-2, tag=size-1)
				reqrecv0.Wait()
				reqsend0.Wait()
				V[0:olap*gsy] = buf
			else:
				reqsend1 = comm.Isend(V[(gsx-2*olap)*gsy:(gsx-olap)*gsy], dest=myrank+1, tag=myrank)
				reqrecv1 = comm.Irecv(buf, source=myrank+1, tag=myrank+1)
				reqrecv1.Wait()
				reqsend1.Wait()
				V[(gsx-olap)*gsy:(gsx)*gsy] = buf

				reqrecv0 = comm.Irecv(buf, source=myrank-1, tag=myrank-1)
				reqsend0 = comm.Isend(V[olap*gsy:2*olap*gsy], dest=myrank-1, tag=myrank)
				reqrecv0.Wait()
				reqsend0.Wait()
				V[0:olap*gsy] = buf


	del RA, BA
	gc.collect()
	return(V)


def Scatter(Vsor, nsize, comm, size, myrank, gsy, olap):

#	print("Nsize out = ", Nsize)
	dspls = [0]
	counts = [nsize[0]*gsy]
	sum = (nsize[0]-2*olap)*gsy
	for i in range(1,size):
		dspls.append(sum)
		counts.append(nsize[i]*gsy)
		sum += (nsize[i]-2*olap)*gsy

	print("size, displacement, counts = ", nsize, dspls, counts)

	if myrank == 0:
#		print("displacement = ", dspls)
		Vin = np.zeros((nsize[0], gsy), dtype=np.float32)
	else:
		Vin = np.zeros((nsize[myrank], gsy), dtype = np.float32)
		Vsor = None

	comm.Scatterv([Vsor,counts,dspls,MPI.FLOAT], Vin, root=0) 

	del Vsor
	gc.collect()

	return Vin

def Gather(Vsor, nsize, comm, size, myrank, gsx, gsy, olap):

	dspls = []
	counts = []
	sum = 0
	for i in range(size):
		dspls.append(sum)
		if i == 0 or i == size -1:
			counts.append((nsize[i]-olap)*gsy)
			sum += counts[-1]
		else:
			counts.append((nsize[i]-2*olap)*gsy)
			sum += counts[-1]

	if myrank == 0:
		Vc = np.zeros((gsx*gsy, 1), dtype=np.float32)
		sendbuf = [Vsor[0:nsize[0]-olap, :], counts[0]]
	elif myrank == size-1:
		Vc = None
		sendbuf = [Vsor[olap:nsize[size-1], :], counts[size-1]]
	else:
		Vc = None
		sendbuf = [Vsor[olap:nsize[myrank]-olap, :], counts[myrank]]

	recvbuf = [Vc,counts,dspls,MPI.FLOAT]
	comm.Gatherv(sendbuf,recvbuf,root=0)

	if myrank == 0 and size != 1:
		Vc = Vc.reshape(gsx, gsy)
#		print(Vc)
	return Vc

Den = []

args = sys.argv
t = 0

print("Initializing")
comm = MPI.COMM_WORLD
size = comm.Get_size()
myrank = comm.Get_rank()
print("My rank = ", myrank)

if myrank == 0:
	print("Nx：" + args[1])
	print("Ny：" + args[2])
	print("Niter_S：" + args[3])
	print("Niter_E：" + args[4])
	print("Maxacc：" + args[5])
	print("Gstartx：" + args[6])
	print("Gbase：" + args[7])
	print("Gth：" + args[8])
	print("Overlap：" + args[9])
	print("Diter：" + args[10])
	print("Plot：" + args[11])
	print("EleFile Name:"+args[12])
	print("RFV:"+args[13])

Nx = int(args[1]); Ny = int(args[2]); Niter = int(args[3]); NiterF = int(args[4]); Maxacc = float(args[5]); 
Gstartx = int(args[6]); Gbase = int(args[7]); Gth = int(args[8]);
Olap = int(args[9]); Diter = int(args[10]); PLT = int(args[11]); Filename = args[12]; RFV = float(args[13]);

#if myrank == 0:
#	InitElec(Filename, Den, RFV)
#	print("Den = ", Den)
InitElec(Filename, Den, Nx, RFV)
	
alldelta = np.zeros(1, np.float32)
teldelta = np.zeros(1, np.float32)

Gsx = Gstartx
Gsy = int(Ny*Gstartx/Nx)
acc = Maxacc

start_time = time.process_time()
if myrank == 0:
	Elec = np.zeros((Gsx, Gsy), np.int8);
	Vin = np.zeros((Gsx, Gsy), np.float32);

	delta = 5
	Electrode(Vin, Elec, Den, Gsx, Gsy, 0, Gsx, t)
	if PLT != 0:
		plotPot(Elec, Gsx, Gsy, 1)                          

	V = np.copy(Vin); Vd = np.copy(Vin);
	while Gsx < Gth:
		print("Gsx, Gsy = ", Gsx, Gsy)

		Vin = Vin.reshape(Gsx*Gsy, 1) 
		V = V.reshape(Gsx*Gsy, 1)    
		Vd = Vd.reshape(Gsx*Gsy, 1)

		RA = BA = []
		RA, BA = SetCSC(Gsx, Gsy, 0)

		Vin = V.copy()
		for iter in range(0, Niter):
			if myrank == 0 and iter%50 < Diter:
				print(iter)
				print("Delta, acc = ", delta, acc)

			# The actual iteration
			V = RA.dot(V)

			Vd = V - Vin
			Rdelta = np.max(abs(Vd))

			V = Vin + Vd*(1+acc)

			Vin = V.copy()
			V = BA.dot(V)
			Vd = V - Vin
			Bdelta = np.max(abs(Vd))

			V = Vin + Vd*(1+acc)
			Vin = V.copy()

			delta = Rdelta
			if delta<Bdelta:
				delta = Bdelta

		Vin = Vin.reshape(Gsx, Gsy)
		del RA, BA, V, Vd, Elec
		gc.collect()
		V = np.zeros((Gsx*Gbase, Gsy*Gbase), np.float32); Vd = np.copy(V);
		for i in range(Gsx):
			for j in range(Gsy):
				V[i*Gbase:(i+1)*Gbase, j*Gbase:(j+1)*Gbase] = Vin[i,j]

		Gsx *= Gbase
		Gsy *= Gbase

		Elec = np.zeros((Gsx, Gsy), np.int8);
		Electrode(V, Elec, Den, Gsx, Gsy, 0, Gsx, t)


		del Vin
		gc.collect()

		Vin = np.copy(V);
		Vd = np.copy(V);

	if PLT != 0:
		plotPot(Vin, Gsx, Gsy, 1)                          
else:
	while Gsx < Gth:
		Gsx *= Gbase
		Gsy *= Gbase

#exit()

Nsize = []
SplitArray(Nsize, size, Gsx, Olap)
Gsxxx = Gsx

if size != 1:
	if myrank != 0:
		V = None
	Vin = Scatter(V, Nsize, comm, size, myrank, Gsy, Olap)

	V = np.copy(Vin); Vd = np.copy(Vin);
	if PLT != 0:
		plotPot(Vin, Nsize[myrank], Gsy, 1)                          

	Elec = np.zeros((Nsize[myrank], Gsy), np.int8);

	offsetx = 0
	for i in range(0, myrank):
		offsetx += Nsize[i]-2*Olap

	Electrode(Vin, Elec, Den, Gsxxx, Gsy, offsetx, Nsize[myrank], t)
	if PLT != 0:
		plotPot(Vin, Nsize[myrank], Gsy, 1)                          

while True:
#	print ("Shape = ", Vin.shape)
	Rdelta, Bdelta, delta, acc = 5, 5, 5, Maxacc

	V = np.copy(Vin);

	V = V.reshape(Nsize[myrank]*Gsy, 1)    

	"""
	Create sparse matrix for SOR method
	"""
	Gsx = Nsize[myrank]
	if size == 1:
		type = 0
	elif myrank == 0:
		type = 1
	elif myrank == size-1:
		type = 2
	else:
		type = 3
		

	RA, BA = SetCSC(Gsx, Gsy, type)
	V = Iteration(V, RA, BA, Nsize, comm, size, myrank, Gsy, Niter, Diter, Olap)

	Vin = np.copy(V)
#	print("V shape = ", Vin.shape)
	Vin = Vin.reshape(Nsize[myrank], Gsy)
	if PLT != 0:
		plotPot(Vin, Gsx, Gsy, 5)
	Nsize_org = Nsize[:]
	print("Nsize_org = ", Nsize_org)

	if Gsxxx*Gbase < Nx:
		Gsxxx *= Gbase
		offsetx_org = 0
		for i in range(0, myrank):
			offsetx_org += Nsize[i]-2*Olap
		V = np.zeros((Nsize_org[myrank]*Gbase, Gsy*Gbase), np.float32)
		SplitArray(Nsize, size, Gsxxx, Olap)

		offsetx = 0
		for i in range(0, myrank):
			offsetx += Nsize[i]-2*Olap

		for i in range(Nsize_org[myrank]):
			for j in range(Gsy):
				V[i*Gbase:(i+1)*Gbase, j*Gbase:(j+1)*Gbase] = Vin[i,j]

		Vin = np.empty((Nsize[myrank], Gsy*Gbase), np.float32)
		Dif = offsetx - offsetx_org*Gbase
		print("offset, offset_org, Dif = ", offsetx, offsetx_org, Dif)
#		print("Vin shape = ", Vin.shape)
		if Dif < 0:
			Dif = 0
		elif Dif+Nsize[myrank] > Nsize_org[myrank]*Gbase:
			Dif = Nsize_org[myrank]*Gbase - Nsize[myrank]
		Vin = V[Dif:Dif+Nsize[myrank] , :]

		V = np.copy(Vin)
		Vd = np.copy(Vin)


		if myrank == 0:
			print("Not Exceed! Gsx, Nx = ", Gsx*Gbase, Nx)
		if PLT != 0:
			plotPot(Vin, Gsx, Gsy, 5)

		Gsx *= Gbase
		Gsy *= Gbase

		Elec = np.zeros((Nsize[myrank], Gsy), np.int8)

		Electrode(Vin, Elec, Den, Gsx, Gsy, offsetx, Nsize[myrank], t)
	

	else:
		if myrank == 0:
			print("Exceed! Gsxxx, Nx = ", Gsxxx*Gbase, Nx)
		if PLT != 0:
			plotPot(Vin, Nsize[myrank], Gsy, 20)

		break;


print
Vc = None
if size != 1:
	Vc = Gather(Vin, Nsize_org, comm, size, myrank, Gsxxx, Gsy, Olap)
else:
	Vc = Vin.reshape(Gsxxx, Gsy)
	


SplitArray(Nsize, size, Nx, Olap)

offsetx_org = 0
for i in range(0, myrank):
	offsetx_org += Nsize_org[i]-2*Olap

offsetx = 0
for i in range(0, myrank):
	offsetx += Nsize[i]-2*Olap


ratioy = Ny/Gsy
ratiox = Nx/Gsxxx

print("offsetx_org, offset, ratiox, ratioy, Nx, Gsxxx", offsetx_org, offsetx, ratiox, ratioy, Nx, Gsxxx)

if myrank == 0:
	V = np.empty((Nx, Ny), np.float32)

	for i in range(Gsxxx-1):
		for j in range(Gsy-1):
			V[int(i*ratiox):int((i+1)*ratiox), int(j*ratioy):int((j+1)*ratioy)] = Vc[i,j]
		V[int(i*ratiox):int((i+1)*ratiox), int((Gsy-1)*ratioy):] = Vc[i,(Gsy-1)]
	for j in range(Gsy-1):
		V[int((Gsxxx-1)*ratiox):, int(j*ratioy):int((j+1)*ratioy)] = Vc[Gsxxx-1,j]

	V[int((Gsxxx-1)*ratiox):, int((Gsy-1)*ratioy):] = Vc[Gsxxx-1,(Gsy-1)]

	del Vc
	gc.collect()

if size != 1:
	if myrank != 0:
		V = None

	Vin = Scatter(V, Nsize, comm, size, myrank, Ny, Olap)

else:
	Vin = np.copy(V)

Vd = np.copy(Vin)

Gsy = Ny

Elec = np.zeros((Nsize[myrank], Gsy), np.int8)
Electrode(Vin, Elec, Den, Nx, Ny, offsetx, Nsize[myrank], t)

V = np.copy(Vin)
Vin = Vin.reshape(Nsize[myrank]*Gsy, 1) 
V = V.reshape(Nsize[myrank]*Gsy, 1)    
Vd = Vd.reshape(Nsize[myrank]*Gsy, 1)


"""
Create sparse matrix for SOR method
"""
if size == 1:
	type = 0
elif myrank == 0:
	type = 1
elif myrank == size-1:
	type = 2
else:
	type = 3

RA, BA = SetCSC(Nsize[myrank], Ny, type)

buf = np.empty((Ny*Olap, 1), np.float32)
Vin = V.copy()

for iter in range(0, NiterF, Diter):
	if myrank == 0:
		if iter%5 < Diter:
			print(iter)
			print("Delta, acc = ", delta, acc)

	# The actual iteration
	for initer in range(Diter):
		V = RA.dot(V)

		Vd = V - Vin
		Rdelta = np.max(abs(Vd))

		V = Vin + Vd*(1+acc)

		Vin = V.copy()
		V = BA.dot(V)
		Vd = V - Vin
		Bdelta = np.max(abs(Vd))

		V = Vin + Vd*(1+acc)
		Vin = V.copy()

#		acc = Maxacc

	'''
	if iter < 1:
		acc = 0.4
	elif iter < 10:
		acc = 0.67
	else:
		Hist = 0.7; Fmax = 0.9
		acc = acc*Hist + (1-Hist)*(Fmax - (Fmax-0.4)*(1/(1+delta/(2*0.005))))
	'''

	if Rdelta > Bdelta:
		delta = Rdelta
	else:
		delta = Bdelta
	teldelta[0] = delta

	comm.Reduce(teldelta, alldelta, op=MPI.MAX, root=0)
	delta = alldelta[0]


	if size != 1:
		if myrank == 0:
			reqsend0 = comm.Isend(V[(Nsize[0]-2*Olap)*Ny:(Nsize[0]-Olap)*Ny], dest=1, tag=0)
			reqrecv0 = comm.Irecv(buf, source=1, tag=1)
			reqsend0.Wait()
			reqrecv0.Wait()
			V[(Nsize[0]-Olap)*Ny:(Nsize[0])*Ny] = buf 
		elif myrank == size-1:
			reqrecv0 = comm.Irecv(buf, source=size-2, tag=size-2)
			reqsend0 = comm.Isend(V[Olap*Ny:(2*Olap)*Ny], dest=size-2, tag=size-1)
			reqrecv0.Wait()
			reqsend0.Wait()
			V[0:Olap*Ny] = buf
		else:
			reqsend1 = comm.Isend(V[(Nsize[myrank]-2*Olap)*Ny:(Nsize[myrank]-Olap)*Ny], dest=myrank+1, tag=myrank)
			reqrecv1 = comm.Irecv(buf, source=myrank+1, tag=myrank+1)
			reqrecv1.Wait()
			reqsend1.Wait()
			V[(Nsize[myrank]-Olap)*Ny:(Nsize[myrank])*Ny] = buf

			reqrecv0 = comm.Irecv(buf, source=myrank-1, tag=myrank-1)
			reqsend0 = comm.Isend(V[Olap*Ny:2*Olap*Ny], dest=myrank-1, tag=myrank)
			reqrecv0.Wait()
			reqsend0.Wait()
			V[0:Olap*Ny] = buf


V = V.reshape(Nsize[myrank], Ny)

Vc = None
if size != 1:
	Vc = Gather(V, Nsize, comm, size, myrank, Nx, Ny, Olap)
else:
	Vc = V

end_time = time.process_time() 
run_time = end_time - start_time	
if myrank == 0:
	print("実行時間：\t", run_time) # 3.1853431999999997[sec]

#	PLT = 1 #Plot V
	if PLT != 0:
		plotPot(Vc, Nx, Ny, 10)
		np.save('field.npy', Vc)
	else:
		np.save('field.npy', Vc)
