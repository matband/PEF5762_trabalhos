import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import os
import numpy as np 


data = pd.read_csv('config.txt', header = None).to_numpy()
data = data[:,0]
Ly = data[0]
Lz = data[1]
ny = int(data[2])
nz = int(data[3])



Y, Z = np.linspace(0, Ly, num=ny), np.linspace(0, Lz, num=nz)
Y, Z = np.meshgrid(Y, Z)
#######################

# Shear Stress
df2 = pd.read_csv(absolute_path+"/shear.csv", sep='\s+',header=None)
df2 = df2.iloc[1:-1,1:-1]
plt.clf()
plt.imshow(df2.to_numpy(), cmap ='Spectral_r' ,interpolation='quadric',extent =[0,Ly,0,Lz])
plt.colorbar(shrink = 0.5)
plt.title("Shear stress colormap at cross section")
plt.show()


#####################
# Warp
df3 = pd.read_csv(absolute_path+"/warp.csv", sep='\s+',header=None)
df3 = df3.iloc[1:-1,1:-1]
plt.clf()
plt.imshow(df3.to_numpy(), cmap ='Spectral_r' ,interpolation='quadric',extent =[0,Ly,0,Lz])
plt.colorbar()
plt.title("Warp function colormap at cross section")
plt.show()
# 3D graph
fig = plt.figure()
ax = plt.axes(projection='3d')
surf = ax.plot_surface( Y[1:-1,1:-1], Z[1:-1,1:-1], df3.to_numpy().T, cmap=cm.Spectral_r,
                       linewidth=0, antialiased=False)
plt.title("Warp surface")
fig.colorbar(surf, shrink = 0.7, location='left')
plt.show()

######################
# Prandtl function
# Colormap
df3 = pd.read_csv(absolute_path+"/prandtl.csv", sep='\s+',header=None)
plt.clf()
plt.imshow(df3.to_numpy(), cmap ='Spectral_r' ,interpolation='quadric',extent =[0,Ly,0,Lz])
plt.colorbar()
plt.title("Prandtl function colormap at cross section")
plt.show()
# 3D graph
fig = plt.figure()
ax = plt.axes(projection='3d')
surf = ax.plot_surface( Y, Z, df3.to_numpy().T, cmap=cm.Spectral_r,
                       linewidth=0, antialiased=False)
plt.title("Membrane analogy at cross section")
fig.colorbar(surf, shrink = 0.7, location='left')
plt.show()

#Convergence
df = pd.read_csv(absolute_path+"/convergence.csv", sep='\s+')
nyaxis= df.columns.values[1:].tolist()
nzaxis = df.values[:, 0]
values = df.values[:, 1:]

# Plot Convergence

# Normalize to [0,1]
norm = plt.Normalize(values.min(), values.max())
colors = cm.cool(norm(values))
rcount, ccount, _ = colors.shape
Y, Z = [float(i) for i in nyaxis],nzaxis
Y, Z = np.meshgrid(Y, Z)

fig = plt.figure()
ax = plt.axes(projection='3d')

surf = ax.plot_surface(Y, Z, values, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
surf.set_facecolor((0,0,0,0))

plt.title("Convergence of Ip for ny x nz nodes")
sct = ax.scatter3D(Y, Z, values, c=values, norm = norm, cmap = 'cool')
fig.colorbar(sct, ax = ax, shrink = 0.5, location='left')
ax.set_xlabel('y nodes')
ax.set_ylabel('z nodes')
plt.show()