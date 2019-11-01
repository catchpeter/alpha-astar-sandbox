"""
check alpha range through thin layers of metal
2019.05.02 pfs

"""


import numpy as np
import pylab as pl

# read in ASTAR txt file
if 0:
	wha = 'gold'
	rho = 19.32 # g/cm3
if 1:
	wha = 'aluminum'
	rho = 2.7 # g/cm3
if 0:
	wha = 'kapton'
	rho = 1.42 # g/cm3
if 0:
	wha = 'liquid-xenon'
	rho = 2.9 # g/cm3
if 1:
	wha = 'silver'
	rho = 10.49 # g/cm3

f = open('alpha-in-%s-electronic-nuclear-total.txt'%wha, 'r')
lines = f.readlines()
m = [] # MeV
ye = [] # MeV cm2 / g
yn = [] # MeV cm2 / g
yt = [] # MeV cm2 / g
for i in range(0,120):
	try:
		m.append( float( lines[i+17][0:9] ) )
		ye.append( float( lines[i+17][10:19] ) )
		yn.append( float( lines[i+17][20:29] ) )
		yt.append( float( lines[i+17][30:39] ) )
	except:
		break
f.close()	
ea = np.asarray(m)*1e3 
se = np.asarray(ye)
sn = np.asarray(yn)
st = np.asarray(yt)
sc = 1e3/1e4 # scale MeV/cm => keV/um

S = st*rho*sc # (MeV cm2/g) * (g/cm3) * (scale) => (keV/um)

pl.figure(8);pl.clf()
pl.plot(ea,se*rho*sc,'bo-',markersize=3,markerfacecolor='w')
pl.plot(ea,sn*rho*sc,'ro-',markersize=3,markerfacecolor='w')
pl.plot(ea,st*rho*sc,'ko-',markersize=3,markerfacecolor='w')
pl.xscale('log')
pl.yscale('log')
pl.axis([1,8000,1,1e4])
pl.title('alpha through %s'%wha)
pl.xlabel('alpha energy / keV')
pl.ylabel('alpha energy loss keV/um')
pl.show()

# E-Z catalog states 100 ug/cm2 gold window over alpha particle sources
# and claims this results in a 15 keV decrease in alpha energy
# 100 ug/cm2 gold => 51.7 nm thickness
# ASTAR implies 22.5 keV loss through this thickness
#
# other E-Z windows:
# 0.9 mg/cm2 aluminized mylar: ~6 microns (not sure of density)
# ah - E-Z seem to claim 6.4 um in the beta particle section


# now I want to integrate the stopping power to see how much energy an 
# alpha has left after traversing a thickness t
# start with a 5.5 MeV alpha

#ea = Ea[Ea<6000]
#st = st[Ea<6000]*rho*sc # stopping power in kev per micron
dx = 0.1 # um
x = np.arange(dx,1000,dx)

# starting conditions:
indx = (np.where([(ea>5000) & (ea<6000)])[1][0])
print('alpha Ei = %d'%ea[indx])
thisE = []; thisE.append(ea[indx]) # 5500 keV
thisS = []; thisS.append(S[indx]) # keV/um
thisT = []; thisT.append(0) # thickness

for step in range(0,5000): # step in loop-land, not physical thickness dx
	if thisE[step]<0:
		break
	kev_lost = thisS[step] * dx #
	thisE.append(thisE[step] - kev_lost) # actual energy loss this step
	thisS.append( np.interp(thisE[step],ea[0:indx],S[0:indx]) ) # interpolated energy loss curve keV/um
	thisT.append(dx+step*dx)
	print('step %d foil thickness %2.1f um alpha %4.1f keV'%(step,dx*step,thisE[step]))
print('need %1.1f um to kill the alpha in %s'%(step*dx,wha))
# shows I need about 10 um to kill a 5.5 MeV alpha in gold
# or about 19 mg/cm2

if 1:
	pl.figure(9);pl.clf()
	pl.plot(thisT,thisE,'k-',label='alpha energy / keV')
	pl.plot(thisT,thisS,'k--',label='keV/um')
	pl.title('alpha through %s'%wha)
	pl.xlabel('thickness traversed (dx=%1.2f um) / um'%dx)
	pl.ylabel('see legend')
	pl.minorticks_on()
	pl.legend()
	pl.show()

