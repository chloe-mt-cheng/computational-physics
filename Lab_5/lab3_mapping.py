#Importing required modules
import struct
import numpy as np
import matplotlib.pyplot as plt
import math

#PROBLEM 1

#PART B
#Opening the data
f = open('N19W156.hgt', 'rb')

#Dummy array of 1442401 zeros
data = np.zeros(1201*1201)

#Using a loop to read the data into an array as done in the lab manual
for i in range(len(data)):
    buffer = f.read(2)
    data[i] = struct.unpack('>h', buffer)[0]

#Reshaping the data to be 1201 by 1201
W = data.reshape(1201,1201)

#Defining ranges of longitude and latitude in degrees
lat = np.linspace(19, 20, 1201)
long = np.linspace(155, 156, 1201)

#Defining constants used in numerical differentiation and finding I
h = 420       #[m]
phi = np.pi   #[rad]
dx = 1        #Diff spacing for x
dy = 1        #Diff spacing for y

#Setting up a 1201x1201 array for dW/dx
dWdx = np.zeros_like(W)

#Using a double for loop to compute dW/dx numerically
for i in range(1201):
    
    #Using forward/backward difference for edge cases
    dWdx[i,0] = (W[i,1] - W[i,0])/(dx)
    dWdx[i,1200] = (W[i,1200] - W[i,1199])/(dx)
    
    #Using central differnece for the rest of the cases
    for j in range(1,1200):
        dWdx[i,j] = (W[i,j+1] - W[i,j-1])/(2*dx)
        
#Setting up a 1201x1201 array for dW/dy
dWdy = np.zeros_like(W)

#Using a double for loop to compute dW/dy numerically
for j in range(1201):
    
    #Using forward/backward difference for edge cases
    dWdy[0,j] = (W[1,j] - W[0,j])/(dy)
    dWdy[1200,j] = (W[1200,j] - W[1199,j])/(dy)
    
    #Using central differnece for the rest of the cases
    for i in range(1,1200):
        dWdy[i,j] = (W[i+1,j] - W[i-1,j])/(2*dy)

#Computing I as in equation 3 and 4 from the lab manual
I = -((np.cos(phi)*dWdx) + (np.sin(phi)*dWdy))/(((dWdx**2) + (dWdy**2) + 1)**(1/2))

#Plotting W as a function of longitude and latitude as a colour map
plt.figure(figsize=(12,10))
plt.imshow(W, extent=[np.amax(long), np.amin(long), np.amin(lat), np.amax(lat)], vmin=0, vmax=4205, aspect='auto', cmap='gist_gray')
plt.title('NASA SRTM Relief Map of Hawaii Island', fontsize=25)
plt.xlabel('Longitude ($x$) [$^{o}$]', fontsize=20)
plt.ylabel('Latitude ($y$) [$^{o}$]', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.colorbar().set_label('Height ($W$) [m]', rotation=90, size=20)
plt.savefig('lab3_q1_plot1.pdf', bbox_inches='tight')
plt.show()

#Plotting I as a function of longitude and latitude as a colour map
plt.figure(figsize=(12,10))
plt.imshow(I, extent=[np.amax(long), np.amin(long), np.amin(lat), np.amax(lat)], vmin=-1, vmax=1, aspect='auto', cmap='gist_gray')
plt.plot(155.590, 19.470, 'rs', markersize=70, fillstyle='none', label='Mauna Loa')
plt.plot(155.458, 19.820, 'bs', markersize=70, fillstyle='none', label='Mauna Kea')
plt.plot(155.840, 19.680, 'gs', markersize=70, fillstyle='none', label='Hualalai')
plt.plot(155.280, 19.410, 'ys', markersize=70, fillstyle='none', label='Kilauea')
plt.title('NASA SRTM Illuminated Relief Map of Hawaii Island', fontsize=25)
plt.xlabel('Longitude ($x$) [$^{o}$]', fontsize=20)
plt.ylabel('Latitude ($y$) [$^{o}$]', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(markerscale=0.22, loc='lower right', fontsize=15)
plt.colorbar().set_label('Surface Illumination ($I$)', rotation=90, size=20)
plt.savefig('lab3_q1_plot2.pdf', bbox_inches='tight')
plt.show()