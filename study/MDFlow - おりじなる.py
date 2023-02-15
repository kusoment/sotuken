""" From "COMPUTATIONAL PHYSICS", 3rd Ed, Enlarged Python eTextBook  
    by RH Landau, MJ Paez, and CC Bordeianu
    Copyright Wiley-VCH Verlag GmbH & Co. KGaA, Berlin;  Copyright R Landau,
    Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2015.
    Support by National Science Foundation"""

'''
MPI multiparticle
Cached Field
Mobility Collision Process
'''
# LaplaceLine.py:  Solve Laplace's eqtn, 3D matplot, close shell to quit

from mpi4py import MPI
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.pylab as p
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
import time
import sys
import gc
import pandas as pd
from numpy.random import *
import os,sys
import mdflow


sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', buffering=1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', buffering=1)
sys.stdin = os.fdopen(sys.stdin.fileno(), 'r', buffering=1)

'''

#グラフ描写



fig = plt.figure(figsize=(6, 10))
ax = fig.add_subplot(4,1,1)
ax.plot(array_t, ztrj, color='r', label='Position Z')
ax2 = ax.twinx()
ax2.plot(array_t, vztrj, label='Velocity Z')
ax.set_yticks(np.arange(0, 3.5E-3*5, step=3.5E-3))
ax.grid(True)
ax.legend()
ax2.legend()

ax3 = fig.add_subplot(4,1,2)
ax3.plot(ztrj,vztrj, label='Velocity R')
ax3.set_xticks(np.arange(0, 3.5E-3*6, step=3.5E-3))
ax3.grid(True)

ax4 = fig.add_subplot(4,1,3)
ax4.plot(array_t, rtrj, color='r', label='Position R')
ax5 = ax4.twinx()
ax5.plot(array_t, vrtrj, label='Velocity R')
ax4.set_ylim([-1.5E-3, 1.5E-3])
ax4.set_yticks(np.arange(-1.5E-3, 1.5E-3, step=0.5E-3))
ax4.grid(True)
ax4.legend()
ax5.legend()

ax6 = fig.add_subplot(4,1,4)
ax6.plot(ztrj,rtrj, label='2D Position')
ax6.set_xticks(np.arange(0, 3.5E-3*6, step=3.5E-3))
ax6.set_yticks(np.arange(0, 1.5E-3, step=0.5E-3))
ax6.grid(True)
ax6.legend()

fig.subplots_adjust(left=0.2, right=0.7, bottom=0.2, top=0.8)
'''

def Cache(VV, Ele, z, r, zorg, rorg, CaZ, CaR):

	Z = int(round(z)); R = int(round(r)); HCyc = SpC//2+1
	print ("z, r, zorg, rorg, nz,nr, CaZ, CaR = ", z, r, zorg, rorg, Nz, Nr, CaZ, CaR)
	MZ = CaZ // 2; LZ = CaZ//4; HZ = (CaZ*3)//4;
	MR = CaR // 2; LR = CaR//4; HR = (CaR*3)//4;
	if z - zorg > LZ:# Z out to position side
		if zorg + HZ < Nz:
			Vwork = np.empty([HCyc, LZ, CaR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([LZ, CaR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(VV[:, LZ:, :], Vwork, axis=1)
			reqsend4.Wait()
			Ele = np.append(Ele[LZ:, :], Elework, axis=0)

			return VV, Ele, Z, rorg

		elif zorg + MZ < Nz:
			LZ = Nz - (zorg+MZ);
			Vwork = np.empty([HCyc, LZ, CaR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([LZ, CaR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(VV[:, LZ:, :], Vwork, axis=1)
			reqsend4.Wait()
			Ele = np.append(Ele[LZ:, :], Elework, axis=0)

			return VV, Ele, Nz-MZ, rorg


	elif zorg - z > LZ:# Z out to negative side
		if zorg - HZ > 0:

			Vwork = np.empty([HCyc, LZ, CaR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([LZ, CaR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(Vwork, VV[:, :HZ, :], axis=1)

			reqsend4.Wait()
			Ele = np.append(Elework, Ele[:HZ, :], axis=0)

			return VV, Ele, Z, rorg

		elif zorg - MZ > 0:

			LZ = (zorg-MZ); HZ = CaZ-LZ
			Vwork = np.empty([HCyc, LZ, CaR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([LZ, CaR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(Vwork, VV[:, :HZ, :], axis=1)

			reqsend4.Wait()
			Ele = np.append(Elework, Ele[:HZ, :], axis=0)

			return VV, Ele, MZ, rorg


	elif r - rorg > LR:# R out to position side
		if rorg + HR < Nr:
			Vwork = np.empty([HCyc, CaZ, LR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([CaZ, LR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			V = np.append(VV[:, :, LR:], Vwork, axis=2)

			reqsend4.Wait()
			Ele = np.append(Ele[:, LR:], Elework, axis=1)

			return VV, Ele, zorg, R

		elif rorg + MR < Nr:
			LR = Nr - (rorg+MR);
			Vwork = np.empty([HCyc, CaZ, LR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([CaZ, LR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(VV[:, :, LR:], Vwork, axis=2)

			reqsend4.Wait()
			Ele = np.append(Ele[:, LR:], Elework, axis=1)

			return VV, Ele, zorg, Nr-MR

	elif rorg - r > LR:# R out to negative side
		if rorg - HR > 0:

			Vwork = np.empty([HCyc, CaZ, LR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([CaZ, LR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(Vwork, VV[:, :, :HR], axis=2)

			reqsend4.Wait()
			Ele = np.append(Elework, Ele[:, :HR], axis=1)

			return VV, Ele, zorg, R
	
		elif rorg - MR > 0:
			LR = (rorg-MR); HR = CaR-LR
			Vwork = np.empty([HCyc, CaZ, LR], dtype = np.float32)
			reqsend3 = comm.Irecv(Vwork, source=0, tag=3)
			Elework = np.empty([CaZ, LR], dtype = np.int8)
			reqsend4 = comm.Irecv(Elework, source=0, tag=4)

			reqsend3.Wait()
			VV = np.append(Vworkz, VV[:, :, :HR], axis=2)

			reqsend4.Wait()
			Ele = np.append(Elework, Ele[:, :HR], axis=1)

			return VV, Ele, zorg, MR


#	print ("Z, R = ", z, r)
	return VV, Ele, zorg, rorg

def diffusion(z, r, D, dt):
	theta=np.pi*rand()
	phi=2*np.pi*rand()
#	rr = abs(randn())*np.sqrt(2*dt*D)
	rr = abs(randn())*0.0001
	x = rr*np.sin(theta)*np.cos(phi)+r
	y = rr*np.sin(theta)*np.sin(phi)
	r = np.sqrt(x*x+y*y)
	z += rr*np.cos(theta)
	return(z, r)

def graP(DZ, DR, cyc, SpC, z, r, OrgZ, OrgR): #Partial Derivative for z of Potential 
#    print("n = ", nz, nr)
	if cyc > SpC // 2:
		cyc = SpC - cyc

	nz = int(z)-OrgZ+CaZ//2
	nr = int(r)-OrgR+CaR//2
	dz = z - int(z)
	dr = r - int(r)

	'''
	difzl = (V[cyc, nz,nr]-V[cyc, nz-1,nr]); difzh = (V[cyc, nz+1,nr]-V[cyc, nz,nr]);
	difrl = (V[cyc, nz,nr]-V[cyc, nz,nr-1]); difrh = (V[cyc, nz,nr+1]-V[cyc, nz,nr]);
	difz = (difzh-difzl)*dz+difzl
	difr = (difrh-difrl)*dr+difrl
#    print("vr = ", vrs, vrl)
	'''
	difzl = DZ[cyc, nz-1,nr]; difzh = DZ[cyc, nz, nr]
	difrl = DR[cyc, nz,nr-1]; difrh = DR[cyc, nz, nr]
	difz = (difzh-difzl)*dz+difzl
	difr = (difrh-difrl)*dr+difrl

	return difz*1E+5, difr*1E+5



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

def InitElec(Filename, Den, Nz, Nr):

	df = pd.read_csv(Filename, index_col=0)
#	print("df = ", df)
	for index, row in df.iterrows():
#		print("index, row = ", index, row)
		Den.extend([[index, float(row['Vol']), float(row['X1'])/Nz, float(row['Y1'])/Nz, float(row['X2'])/Nz, float(row['Y2'])/Nz]])
#	print("Den = ", Den)
	del df
	gc.collect()

def InitParticle(Filename, Particle, Nz, Nr):

	df = pd.read_csv(Filename, index_col=0)
#ParticleNumber				0
#Z							1
#R							2
#CaZ						3
#CaR						4
#Step						5
#R/nm						6
#Q							7
#ParticleDensity/g/cm3		8
#ParticleMass/MU			9
#ParticleRelaxationTime/s	10
#ParticleFreeTime/s			11
#ParticleMobility/Vs/m^2	12
#StokesMobility				13
#CorrelatedMobility			14
#Cunningham_A				15
#B							16
#C							17
#SelectMobility/012			18
#GasMass/MU					19
#GasRadius/nm				20
#GasPressure/Pa				21
#GasTemperature/K			22
#GasFreetime/s				23
#Viscosity/Pas				24
#GasFreePath/nm				25
#Gravity					26
#GravityAngle/degree		27
#Flow/m/s					28
#FlowAngle/degree			29

	for index, row in df.iterrows():
#		print("index, row = ", index, row)
		Particle.extend([[index, float(row['Z']), float(row['R']+Nr//2), int(row['CaZ']), int(row['CaR']), int(row['Step']), float(row['R/nm']), float(row['Q']), float(row['ParticleDensity/g/cm3']), float(row['ParticleMass/MU']), float(row['ParticleRelaxationTime/s']), float(row['ParticleFreeTime/s']), float(row['ParticleMobility/Vs/m^2']), float(row['StokesMobility']), float(row['CorrelatedMobility']), float(row['Cunningham_A']), float(row['B']), float(row['C']), int(row['SelectMobility/012']), float(row['GasMass/MU']), float(row['GasRadius/nm']), float(row['GasPressure/Pa']), float(row['GasTemperature/K']), float(row['GasFreeTime/s']), float(row['Viscosity/Pas']), float(row['GasFreePath/nm']), float(row['Gravity']), float(row['GravityAngle/degree']), float(row['Flow/m/s']), float(row['FlowAngle/degree']), int(round(float(row['Z']))), int(round(float(row['R']+Nr//2)))]])
	print("Particle = ", Particle)
	del df
	gc.collect()

def Electrode(Elec, Den, Nz, Nr, offsetx, deltax):

#	print(Den)
	for e in Den:
		i, V, x1, y1, x2, y2  = e
#		if myrank == 0:
#			print("e, x1, y1 ... = ", e, x1, y1, x2, y2, i, V)
			
		x1 = int(x1)	#Always lower left
		y1 = int(y1)
		x2 = int(x2)	#Always upper right
		y2 = int(y2)
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
			elif 0 <= y2 < Nr//2+1:
				y1 = 0
				y2 = y2 + dy
			else:
				y1 = 0
				y2 = Nr//2+1
		elif 0 <= y1 < Nr//2+1:
			if 0 <= y2 < Nr//2:
				y2 = y2 + dy
			else:
				y2 = Nr//2
		else:
			continue

		Elec[x1:x2, Nr//2-y1:Nr//2-y2] = i #Lower Half Electrode Number is i < 256 Center R = Nr//2
		Elec[x1:x2, y1+Nr//2:y2+Nr//2] = i #Upper Half Electrode Number is i < 256
#		Pot[x1:x2, y1:y2] = V

#		if myrank == 0:
#			print("Elec = ", Elec)
#			print("Pot = ", Pot)

def Iteration(myrank, size, comm, trj, Maxiter, Hz, SpC, ChS, SvS):

	#VZ += 10/3.5E-3

	vz = trj.iloc[0]['vz'] #軸初期速度
	vr = trj.iloc[0]['vr'] #半径初期速度
	z = trj.iloc[0]['z'] #軸初期位置
	r = trj.iloc[0]['r']+Nr//2 #半径初期位置
	t = trj.iloc[0]['t'] #初期時間
	zo_save = z; ro_save = r # 記録位置
#	h = 0.0005 #ルンゲクッタ法の刻み幅
	h = 1/(Hz*SpC) #ルンゲクッタ法の刻み幅
	HCyc = SpC//2+1

	Pend = 29 #Particle Parameter End

#	K=3.36E-6
	#ルンゲクッタ法のメイン計算

	
	pos = np.empty(2, dtype = np.float32)
	end = np.empty(1, dtype = np.int8)

	if myrank == 0:
		for i in range(1, size):
			Z = Particle[i-1][Pend+1]; R = Particle[i-1][Pend+2]; #Particle Current Position in Array
			CaZ = Particle[i-1][3]; CaR = Particle[i-1][4];
			MZ = CaZ//2; MR = CaR//2
			print("i, center, cazr, Nzr = ", i, Z, R, CaZ, CaR, Nz, Nr)

			LZ = Z - MZ; HZ = Z + MZ; LR = R-MR; HR = R+MR
			if Z-MZ < 0:
				LZ = 0; HZ = CaZ; Z = MZ
			elif Z + MZ > Nz -1:
				LZ = Nz - 1 - CaZ; HZ = Nz - 1; Z = HZ-MZ
			if R-MR < 0:
				LR = 0; HR = CaR; R = MR
			elif R + MR > Nr -1:
				LR = Nr - 1 - CaR; HR = Nr - 1; R = HR - MR

			print("Z, R, LZ, HZ, LR, HR = ", Z, R, LZ, HZ, LR, HR)
			Particle[i-1][Pend+1]=Z; Particle[i-1][Pend+2]=R; #Particle Current Position in Array

			Vwork = V[:, LZ:HZ, LR:HR]
			Vwork = Vwork.reshape(1,HCyc*CaZ*CaR)
			reqsend0 = comm.Isend(Vwork, dest=i, tag=0)


			Elework = np.empty([CaZ, CaR], dtype = np.int8)
			Elework = Elec[LZ:HZ, LR:HR]
			Elework = Elework.reshape(CaZ*CaR)
			reqsend2 = comm.Isend(Elework, dest=i, tag=2)

			pos[0] = Z; pos[1] = R
			reqsend10 = comm.Isend(pos, dest=i, tag=10)
			reqsend0.Wait()
			reqsend2.Wait()
			reqsend10.Wait()

		reqpos = []; reqend = []
		alive = list(range(1, size))
		for i in range(1, size):
#			print("Requesting Cache and End Search for rank = ", i)
			reqpos.append(comm.Irecv(pos, source=i, tag=1000))
			reqend.append(comm.Irecv(end, source=i, tag=2000))

#		print("Finish Search & Alive =!", alive)
		while(len(alive) != 0):
			for i in alive:
				if reqpos[i-1].Test():
					reqpos[i-1].Wait()
#					print ("rank, pos = ", i, pos)
					z = pos[0]; r = pos[1]
					Particle[i-1][1]=z; Particle[i-1][2]=r
					Z = int(round(z)); R = int(round(r))
					OrgZ = Particle[i-1][Pend+1]; OrgR = Particle[i-1][Pend+2]; #Cache Origin Revised
					CaZ = Particle[i-1][3]; CaR = Particle[i-1][4]
					MZ = CaZ//2; MR = CaR//2
					HZ = (CaZ * 3) // 4; LZ = CaZ // 4
					HR = (CaR * 3) // 4; LR = CaR // 4

					if z - OrgZ >= LZ:# Z out to positive side
						if OrgZ + HZ < Nz:
							print("Source Z out to positive side ID, pos, org, Cachesize, //4= ", i-1, pos, OrgZ, OrgR, CaZ, CaR, LZ, LR)
							Vwork = V[:, OrgZ+MZ:OrgZ+MZ+LZ, OrgR-MR:OrgR+MR]
#							print("Vwork.shape and Numbers= ", Vwork.shape, (HCyc)*(LZ+1)*(CaR+1))
							Vwork = Vwork.reshape(1,HCyc*LZ*CaR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ+MZ:OrgZ+MZ+LZ, OrgR-MR:OrgR+MR]
							Elework = Elework.reshape(1,LZ*CaR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+1]=Z;
#							print("Finish Sending Z out positive & position = ", pos)
						elif OrgZ + MZ < Nz:
				#			exit()
							LZ = Nz-1 - (OrgZ+MZ);
							Vwork = V[:, OrgZ+MZ:OrgZ+MZ+LZ, OrgR-MR:OrgR+MR]
							Vwork = Vwork.reshape(1,HCyc*LZ*CaR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ+MZ:OrgZ+MZ+LZ, OrgR-MR:OrgR+MR]
							Elework = Elework.reshape(1,LZ*CaR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+1]=Nz-MZ;
#							print("Finish Sending Z out positive & position = ", pos)

					elif OrgZ - z >= LZ:# Z out to negative side
						if OrgZ - HZ > 0:
							print("Source Z out to negative side ID, pos, org, Cachesize, //4 = ", i-1, pos, OrgZ, OrgR, CaZ, CaR, LZ, LR)
							Vwork = V[:, OrgZ-MZ-LZ:OrgZ-MZ, OrgR-MR:OrgR+MR]
							Vwork = Vwork.reshape(1,HCyc*LZ*CaR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ-LZ:OrgZ-MZ, OrgR-MR:OrgR+MR]
							Elework = Elework.reshape(1,LZ*CaR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+1]=Z;
#							print("Finish Sending Z out negative & position = ", pos)

						elif OrgZ - MZ > 0:
							LZ = (OrgZ-MZ); HZ = CaZ-LZ
							Vwork = V[:, OrgZ-MZ-LZ:OrgZ-MZ, OrgR-MR:OrgR+MR]
							Vwork = Vwork.reshape(1,HCyc*LZ*CaR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ-LZ:OrgZ-MZ, OrgR-MR:OrgR+MR]
							Elework = Elework.reshape(1,LZ*CaR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+1]=MZ;

					elif r - OrgR >= LR:# R out to position side
						if OrgR + HR < Nr:
							print("Source R out to positive side ID, pos, org, Cachesize, //4 = ", i-1, pos, OrgZ, OrgR, CaZ, CaR, LZ, LR)
							Vwork = V[:, OrgZ-MZ:OrgZ+MZ, OrgR+MR:OrgR+MR+LR]
							Vwork = Vwork.reshape(1,HCyc*CaZ*LR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ:OrgZ+MZ, OrgR+MR:OrgR+MR+LR]
							Elework = Elework.reshape(1,CaZ*LR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+2]=R
#							print("Finish Sending R out positive & position = ", pos)
						elif OrgR + MR < Nr:
							LR = Nr- (OrgR+MR)
							Vwork = V[:, OrgZ-MZ:OrgZ+MZ, OrgR+MR:OrgR+MR+LR]
							Vwork = Vwork.reshape(1,HCyc*CaZ*LR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ:OrgZ+MZ, OrgR+MR:OrgR+MR+LR]
							Elework = Elework.reshape(1,CaZ*LR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+2]=Nr-MR

					elif OrgR - r >= LR:# R out to negative side
						if OrgR - HR > 0:
							print("Source R out to negative side ID, pos, org, Cachesize, //4 = ", i-1, pos, OrgZ, OrgR, CaZ, CaR, LZ, LR)
							Vwork = V[:, OrgZ-MZ:OrgZ+MZ, OrgR-MR-LR:OrgR-MR]
							Vwork = Vwork.reshape(1,HCyc*CaZ*LR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ:OrgZ+MZ, OrgR-MR-LR:OrgR-MR]
							Elework = Elework.reshape(1,CaZ*LR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+2]=R
						elif OrgR - MR > 0:
							LR = (OrgR-MR); HR = CaR-LR
							Vwork = V[:, OrgZ-MZ:OrgZ+MZ, OrgR-MR-LR:OrgR-MR]
							Vwork = Vwork.reshape(1,HCyc*CaZ*LR)
							reqsend3 = comm.Isend(Vwork, dest=i, tag=3)
							Elework = Elec[OrgZ-MZ:OrgZ+MZ, OrgR-MR-LR:OrgR-MR]
							Elework = Elework.reshape(1,CaZ*LR)
							reqsend4 = comm.Isend(Elework, dest=i, tag=4)
							reqsend3.Wait()
							reqsend4.Wait()
							Particle[i-1][Pend+2]=MR

					reqpos[i-1] = comm.Irecv(pos, source=i, tag=1000) # next page request
#					print("Master Req position")

#				print("Master End Get = ", reqpos[i-1].Test())
#				print("Master End Get = ", reqend[i-1].Test())

				if reqend[i-1].Test() == True:
#					print ("rank, end = ", i, end)
#					print("Master End Get = ", reqend[i-1].Test())
					reqend[i-1].Wait()
					reqpos[i-1].Cancel()
					alive.remove(i)


	else:
#EleNum						0
#Z							1
#R							2
#CaZ						3
#CaR						4
#Step						5
#R/nm						6
#Q							7
#ParticleDensity/g/cm3		8
#ParticleMass/MU			9
#ParticleRelaxationTime/s	10
#ParticleFreeTime/s			11
#ParticleMobility/Vs/m^2	12
#StokesMobility				13
#CorrelatedMobility			14
#Cunningham_A				15
#B							16
#C							17
#SelectMobility/012			18
#GasMass/MU					19
#GasRadius/nm				20
#GasPressure/Pa				21
#GasTemperature/K			22
#GasFreetime/s				23
#Viscosity/Pas				24
#GasFreePath/nm				25
#Gravity					26
#GravityAngle/degree		27
#Flow/m/s					28
#FlowAngle/degree			29

		CenterZ = float(Particle[myrank-1][1]); CenterR = float(Particle[myrank-1][2]);
		CaZ = Particle[myrank-1][3]; CaR = Particle[myrank-1][4]; Step = Particle[myrank-1][5]
#		Mob = float(Particle[myrank-1][3]); 
		PRad = float(Particle[myrank-1][6]); 
		Q = float(Particle[myrank-1][7]);
		Density = float(Particle[myrank-1][8]); 
		Pmass = 4*np.pi/3*(PRad*1E-9)**3*Density*1000 #in kg
		Gmass = float(Particle[myrank-1][19])*1.6605390666E-27 # in kg 
		GRad = float(Particle[myrank-1][20]);
		Pressure = float(Particle[myrank-1][21]); Temperature = float(Particle[myrank-1][22]);
		OrgZ = Particle[myrank-1][Pend+1]; OrgR = Particle[myrank-1][Pend+2]
		VisSel = Particle[myrank-1][18]

		CunA = float(Particle[myrank-1][15])
		CunB = float(Particle[myrank-1][16])
		CunC = float(Particle[myrank-1][17])

		G = float(Particle[myrank-1][26])
		Gang = float(Particle[myrank-1][27])
		Flow = float(Particle[myrank-1][28])
		Fang = float(Particle[myrank-1][29])


#		D = Mob*1.3806503E-023/(1.60217653E-019*Q)

		Omega = np.pi*((PRad+GRad)*1E-9)**2
		Gdensity = 6.0221415E+023*Pressure/(8.314472*Temperature)
		Gvel = (8*1.3806503E-023*Temperature/(Gmass*np.pi))**0.5;
		tau = 1/(Omega*(8*1.3806503E-023*Temperature/(Gmass*np.pi))**0.5*Gdensity);
		hh = h/Step
		VisMD = 4/3*Gmass/(Pmass+Gmass) /tau *hh
		MobMD = 3/16*Q*1.60217662E-019/Gdensity*(1/Pmass+1/Gmass)**0.5*(2*np.pi/(1.3806503E-023*Temperature))**0.5/Omega#=S2*1.6605390666E-27*(8*1.3806503E-23*$V2/($S2*1.6605390666E-27*PI()))^0.5/(2*SQRT(2)*PI()*(2*T2)^2*0.000000000000000001)
		Gfreetime = 1/np.sqrt(2.0)*1/(4*np.pi*GRad**2*1E-18*Gvel*Gdensity)
#=W2*SQRT(8*1.3806503E-23*V2/(PI()*S2*1.6605390666E-27))*1000000000
		Gfreepath = Gfreetime*Gvel
#=S2*1.6605390666E-27*(8*1.3806503E-23*$V2/($S2*1.6605390666E-27*PI()))^0.5/(2*SQRT(2)*PI()*(2*T2)^2*0.000000000000000001)
		Gvis = Gmass*Gvel/(2*np.sqrt(2.0)*np.pi*(2*GRad*1E-9)**2)
#=H2*1.60217662E-19/(6*PI()*X2*G2*0.000000001)
		MobST = Q*1.60217662E-019/(6*np.pi*Gvis*PRad*1E-9)
#=N2*(1+(0.5*Y2/G2)*(P2+Q2*EXP(R2/(0.5*Y2/G2))))
		Knun = Gfreepath/(2*PRad*1E-9)
		MobCR = MobST*(1+Knun*(CunA+CunB*np.exp(CunC/Knun)))
		VisST = Q*1.60217662E-019/MobST*hh/Pmass
		VisCR = Q*1.60217662E-019/MobCR*hh/Pmass
		if(VisSel == 0):
			Vis = VisMD
		elif (VisSel == 1):
			Vis = VisST
		elif (VisSel == 2):
			Vis = VisCR
#		Vis = 0

		print("tau, Vis = ", tau, Vis)
		print("Q, PRad = ", Q, PRad)
		print("Gft, Gfp, Gvis, MobST, MobCR = ", Gfreetime, Gfreepath, Gvis, MobST, MobCR)
		print("VisMD, VisST, VisCR = ", VisMD, VisST, VisCR)
		CaZH = CaZ//2; CaRH = CaR//2
#		OrgZ, OrgR = round(CenterZ), round(CenterR)
		it = 0; cyc = 0
		svs = 0
		Pcy = 0; 
		Pcyn = 1
#		acc = -Mob*h*1E+5*1E+5 # meter -> grid size of 0.01 mm & Gradient and Direction of Move

		acc = -Q*1.60217662E-019/Pmass*1E+10*hh*hh
		gx = np.cos(Gang/180.0*np.pi)*G; gy = np.sin(Gang/180.0*np.pi)*G;
		gx = gx*1E+5*hh*hh; gy=gy*1E+5*hh*hh
		fx = np.cos(Fang/180.0*np.pi)*Flow; fy = np.sin(Fang/180.0*np.pi)*Flow;
		fx = fx*1E+5*hh; fy=fy*1E+5*hh

		print("acc, Vis, K = ", acc, Vis, MobMD)
		OZ = -OrgZ+CaZH
		OR = -OrgR+CaRH
		print("OZ, OR = ", OZ, OR)
#		print("center, cazr, Nzr = ", CenterZ, CenterR, CaZ, CaR, Nz, Nr)

		VV = np.empty([HCyc, CaZ, CaR], dtype=np.float32)
		Ele = np.empty([CaZ, CaR], dtype=np.int8)

		reqrecv0 = comm.Irecv(VV, source=0, tag=0)
		reqrecv2 = comm.Irecv(Ele, source=0, tag=2)
		reqrecv10 = comm.Irecv(pos, source=0, tag=10)
		reqrecv0.Wait()
		reqrecv2.Wait()
		reqrecv10.Wait()
		OrgZ = int(pos[0]); OrgR = int(pos[1]) #Already Origin Positions are integers
		OZ = -OrgZ+CaZH
		OR = -OrgR+CaRH
		print("OZ, OR = ", OZ, OR)

#		print("DZ, DR, Ele shape = ", DZ.shape, DR.shape, Ele.shape)


		posp = mdflow.Position(z, r, vz, vr, 0, 0, OZ, OR)
		VV = VV.reshape(1,HCyc*CaZ*CaR)

		conp = mdflow.Condition(Pcy, Pcyn, Step, CaZ, CaR, acc, Vis, gx, gy, fx, fy) #conditions
		posp.show()
		conp.show()
#		print("z, r, Pcy, Pcyn = ", posp.z, posp.r, conp.cur, conp.nex)
		for it in range(Maxiter):
#			print("before mdcal\n")
			posp = mdflow.mdcal(posp, conp, VV)
#			posp.show()
			cyc += 1
			t += h
			nz = posp.nz; nr = posp.nr;
#			dvz = difzf; dvr = difrf; tmp_se = pd.Series( [t, z, r, vz, vr, dvz, dvr], index=trj.columns )
#			trj = trj.append(tmp_se,  ignore_index=True )
			if 2*cyc > SpC:
				Pcy = SpC - cyc

				if cyc == SpC:
					cyc = 0
					if Ele[nz, nr] != 0:
						break;
					svs += 1
					if svs % ChS == 0:
						if posp.r>Nr-2:
							posp.r = Nr-2

						if (abs(posp.z - OrgZ) > CaZ//4 and posp.z > CaZ//2 and posp.z < Nz-CaZ//2-1) or (abs(posp.r - OrgR) > CaR//4 and posp.r > CaR//2 and posp.r < Nr-CaR//2-1):
							pos[0] = posp.z; pos[1] = posp.r
							comm.Isend(pos, dest=0, tag=1000) # request potential

							VV = VV.reshape(HCyc, CaZ, CaR) 

							VV, Ele, OrgZ, OrgR = Cache(VV, Ele, posp.z, posp.r, OrgZ, OrgR, CaZ, CaR)


							VV = VV.reshape(1,HCyc*CaZ*CaR)

							posp.OZ = -OrgZ+CaZH
							posp.OR = -OrgR+CaRH

					if svs % SvS == 0:
						dvz = (posp.z - zo_save)*Hz/SvS ; dvr = (posp.r - ro_save)*Hz/SvS;
						zo_save = posp.z; ro_save = posp.r;
						tmp_se = pd.Series( [t, posp.z, posp.r-Nr//2, posp.vz, posp.vr, dvz, dvr], index=trj.columns )
						trj = trj.append(tmp_se,  ignore_index=True )


			else:
				Pcy = cyc
			if 2*(cyc+1) > SpC:
				Pcyn = SpC - (cyc+1)
			else:
				Pcyn = cyc + 1
			conp.cur = Pcy; conp.nex = Pcyn



		return(trj)


args = sys.argv
t = 0

print("Initializing")
comm = MPI.COMM_WORLD
size = comm.Get_size()
myrank = comm.Get_rank()
print("My rank = ", myrank)

if myrank == 0:
	print("Directory Name:"+args[1])
	print("ParticleFile Name:"+args[2])
	print("ElectrodeFile Name:"+args[3])
	print("Nz：" + args[4])
	print("Nr：" + args[5])
	print("Hz:" + args[6])
	print("Step/Cycle:" + args[7])
	print("Niter：" + args[8])
	print("CheckStep：" + args[9])
	print("SaveStep：" + args[10])
	print("Plot?：" + args[11])

DirN = args[1]
Particlename = DirN+'/'+args[2]; Electrodename = DirN+'/'+args[3]; Nz = int(args[4]); Nr = int(args[5])*2-1; 
Hz = float(args[6]); SpC = int(args[7]); Niter = int(args[8]); ChS = int(args[9]); SvS = int(args[10])
PLT = int(args[11])

#V = np.load(args[1])

if myrank == 0:
	name = DirN +'/'+ str(0) + '-' + str(SpC) + '.npy'
	V = np.load(name)
	Vtemp = np.copy(V)
	Vtemp = Vtemp[:, ::-1]
	V = np.append(Vtemp, V[:,0:Nr//2], axis=1)

	for i in range(1, SpC//2+1):
		name = DirN +'/'+ str(i) + '-' + str(SpC) + '.npy'
		Vtemp = np.load(name)
		Vtemp2 = np.copy(Vtemp)
		Vtemp2 = Vtemp2[:, ::-1]
		Vtemp = np.append(Vtemp2, Vtemp[:, 0:Nr//2], axis=1)
		V = np.vstack([V, Vtemp])

	V = V.reshape(SpC // 2+1, Nz, Nr)
	print("Vshape = ", V.shape)


#V = np.load(args[1])
alldelta = np.zeros(1, np.float32)
teldelta = np.zeros(1, np.float32)

'''
allparam = np.zeros(6, np.float32)
telparam = np.zeros(6, np.float32)

if myrank == 0:
	allpos = np.zeros([2*(size-1)], dtype=np.float32)
	telpos = None
else:
	allpos = None
	telpos = [Center[myrank, 0], Center[myrank, 1]]

comm.Gather(telpos,allpos,root=0)

if myrank == 0 and size != 1:
	Vc = Vc.reshape(gsx, gsy)
	print(Vc)
'''

end = np.empty(1, np.int8)
end[0] = 111
'''
sendbuf = None
if rank == 0:
    sendbuf = np.empty([size, 100], dtype='i')
    sendbuf.T[:,:] = range(size)
recvbuf = np.empty(100, dtype='i')
comm.Scatter(sendbuf, recvbuf, root=0)
assert np.allclose(recvbuf, rank)
'''

Den = []
Particle = []

if myrank == 0:
	InitParticle(Particlename, Particle, Nz, Nr)

	InitElec(Electrodename, Den, Nz, Nr)
	Elec = np.zeros((Nz, Nr), np.int8);
	Electrode(Elec, Den, Nz, Nr, 0, Nz)
	trj = pd.DataFrame([[0, 0, 0, 0, 0, 0, 0]], columns = ['t', 'z', 'r', 'vz', 'vr', 'dvz', 'dvr'], dtype = np.float64)
	print("Source tjr = ", trj)
	trj = Iteration(myrank, size, comm, trj, Niter, Hz, SpC, ChS, SvS)
	teldelta[0] = 0
	comm.Reduce(teldelta, alldelta, op=MPI.MAX, root=0)
	delta = alldelta[0]

else:

	InitParticle(Particlename, Particle, Nz, Nr)

	CenterZ = float(Particle[myrank-1][1]); CenterR = float(Particle[myrank-1][2]);


	trj = pd.DataFrame([[0, CenterZ, CenterR-Nr//2, 0, 0, 0, 0]], columns = ['t', 'z', 'r', 'vz', 'vr', 'dvz', 'dvr'], dtype = np.float64)
	print(trj)

	start_time = time.process_time()

#K = 3.36E-5
	trj = Iteration(myrank, size, comm, trj, Niter, Hz, SpC, ChS, SvS)
	end_time = time.process_time() 
	run_time = end_time - start_time
	teldelta[0] = run_time
	reqend = comm.Isend(end, dest=0, tag=2000)
	reqend.Wait()

	comm.Reduce(teldelta, alldelta, op=MPI.MAX, root=0)
	delta = alldelta[0]

if myrank == 0:
	print("実行時間：\t", delta) # 3.1853431999999997[sec]

if myrank != 0:
	print("Trajetory of ID ", myrank)
	print(trj)
	trjdf = pd.DataFrame(trj)
	trjdf.to_csv(str(myrank-1)+".csv")

	if PLT != 0:
		ax = trj.plot(kind='scatter', x='z', y='r')
		plt.show()


'''
fig = plt.figure(figsize=(6, 10))
ax = fig.add_subplot(4,1,1)
ax.plot(trj[t].t, trj[z].t, color='r', label='Position Z')
ax2 = ax.twinx()
ax2.plot(trj['t'].t, trj['vz'].t, label='Velocity Z')
#ax.set_yticks(np.arange(0, 3.5E-3*5, step=3.5E-3))
ax.grid(True)
ax.legend()
ax2.legend()
'''
'''
ax3 = fig.add_subplot(4,1,2)
ax3.plot(ztrj,vztrj, label='Velocity R')
ax3.set_xticks(np.arange(0, 3.5E-3*6, step=3.5E-3))
ax3.grid(True)

ax4 = fig.add_subplot(4,1,3)
ax4.plot(array_t, rtrj, color='r', label='Position R')
ax5 = ax4.twinx()
ax5.plot(array_t, vrtrj, label='Velocity R')
ax4.set_ylim([-1.5E-3, 1.5E-3])
ax4.set_yticks(np.arange(-1.5E-3, 1.5E-3, step=0.5E-3))
ax4.grid(True)
ax4.legend()
ax5.legend()

ax6 = fig.add_subplot(4,1,4)
ax6.plot(ztrj,rtrj, label='2D Position')
ax6.set_xticks(np.arange(0, 3.5E-3*6, step=3.5E-3))
ax6.set_yticks(np.arange(0, 1.5E-3, step=0.5E-3))
ax6.grid(True)
ax6.legend()
'''
#	PLT = 1 #Plot V
#	if PLT != 0:
#		plotPot(Vc, Nx, Ny, 10)
#		np.save('field.npy', Vc)
#	else:
#		np.save('field.npy', Vc)
