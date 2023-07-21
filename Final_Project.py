# -*- coding: utf-8 -*-
"""
Created on Wed May  5 19:38:30 2021

@author: Erik
"""
import Semiamp
import Vmod
import numpy as np
import matplotlib.pyplot as plt  
from astropy.io import ascii

Data=ascii.read('Combined_Data_HARPS_PFS.csv', names=('date','vel','errvel'))
obstimes = Data['date']

# physical and system constants
inc = 71.8 #inclination of orbit
epoch = 2458382.228
G = 6.6738e-11 #Gravitational Constant (Nm^2/kg^2)
M = 1.12 * 1.9891e30 #Mass of the star (kg)
m = 5.9722e24 #Mass of the Earth (kg)
GM = G * M
Gm = G * m
P = 2.476 #period
d = np.cbrt(((P**2) * (G * (M+m))) / (4 * (np.pi**2))) #distance between planet and star
b1 = d / (1 + (M/m)) #distance from star to barycenter
b2 = d - b1 #distance from planet to barycenter
print('d =', d, 'b1 =', b1, 'b2 =', b2)
for i in range(1, 41):
    m = i / 2 #mass of planet in earth masses
    K = Semiamp.Semiamp(M, m, inc, G, P)
    #print(K)
    
    obstimes = Data['date']
    injectedRVs = np.zeros(len(obstimes))
    
    for j in range(0, len(obstimes)):
        injectedRVs[j] = Vmod.Vmod(K, P, epoch, obstimes[j], offset = False) + Data['vel'][j]
        # print(Data['date'][j], injectedRVs[j], Data['errvel'][j])

# acceleration of the star
# def acceleration1(r1):
#     x = r1[0]
#     y = r1[1]
#     a = ((x**2)+(y**2))**0.5
#     return np.array([-(GM)*(x/(a**3)), -(GM)*(y/(a**3))])

def acceleration1(r1):
    
    x=r1[0]        #unpack the array r into local variables
    y=r1[1]
    radius_cubed = (x*x + y*y)**(3/2)
    
    # this returns the acceleration vector
    return -GM*r1/radius_cubed
#acceleration of the planet
# def acceleration2(x, y):
#     x = r2[0]
#     y = r2[1]
#     a = ((x**2)+(y**2))**0.5
#     return np.array([-(Gm)*(x/(a**3)), -(Gm)*(y/(a**3))])



# start = obstimes[0]
# end = obstimes[-1] #number of seconds in two revolutions around the sun
start = 0
end = 63e6
numSteps = 2000
stepSize = (end-start)/numSteps #makes time step one hour
print(stepSize)
#velocities
RV = np.mean(injectedRVs)   #radial velocity of the star
print('RV =', RV)
# vel = np.mean(Data['vel']) #velocity of the planet
# RV = 3.0287e4
# vel = 3.0287e4




tpoints = np.arange(start,end,stepSize)

x1pts = []
y1pts = []
# x1velpts = []
# y1velpts = []

# x2pts = []
# y2pts = []
# x2velpts = []
# y2velpts = []

# initial conditions for acceleration and velocity; summarized in vector r
# r1 = np.array([b1, 0, 0, RV], float)  
# r2 = np.array([b2, 0, 0, vel], float) 

# current_accel2 = acceleration2(r2[0], r2[3])
r1 = np.array([0.0,b1],float)
v1 = np.array([RV,0.0],float)
current_accel1 = acceleration1(r1)

for t in tpoints:
    x1pts.append(r1[0])
    y1pts.append(r1[1])
    # x1velpts.append(r1[2])
    # y1velpts.append(r1[3])
    
    v1 += current_accel1*stepSize/2      # update velocity half-way, acceleration saved from before
    r1 += v1*stepSize
    current_accel1 = acceleration1(r1)    # update acceleration
    v1 += current_accel1*stepSize/2      # update velocity to the end of time step
    
    # r1[2] += current_accel1[0]*stepSize/2 # update x-velocity
    # r1[3] += current_accel1[1]*stepSize/2 # update y-velocity
    # r1[0] += r1[2]*stepSize # update x coordinate using x-velocity
    # r1[1] += r1[3]*stepSize # same for y coordinate
    # current_accel1 = acceleration1(r1[0],r1[1]) # update acceleration using x and y coords
    # r1[2] += current_accel1[0]*stepSize/2 # updating velocities again
    # r1[3] += current_accel1[1]*stepSize/2
    # print('r1 =', r1)
    # print('current accel =', current_accel1)
# for t in tpoints:
#     x2pts.append(r2[0])
#     y2pts.append(r2[1])
#     x2velpts.append(r2[2])
#     y2velpts.append(r2[3])

#     r2[2] += current_accel2[0]*stepSize/2 # update x-velocity
#     r2[3] += current_accel2[1]*stepSize/2 # update y-velocity
#     r2[0] += r2[2]*stepSize # update x coordinate using x-velocity
#     r2[1] += r2[3]*stepSize # same for y coordinate
#     current_accel2 = acceleration2(r2[0],r2[1]) # update acceleration using x and y coords
#     r2[2] += current_accel2[0]*stepSize/2 # updating velocities again
#     r2[3] += current_accel2[1]*stepSize/2


# plot results
plt.plot(x1pts, y1pts, 'b-', label = 'Star Orbiting barycenter')
# plt.plot(x2pts, y2pts, 'g-', label = 'Planet Orbiting barycenter')
plt.plot(0, 0, 'yo', label = 'Barycenter')
plt.title('Earth Orbiting the Sun')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best')
# plt.gca().set_aspect('equal', adjustable = 'box')
plt.show()