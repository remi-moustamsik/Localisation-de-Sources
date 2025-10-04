import numpy as np
import matplotlib.pyplot as plt
from math import *

c = 340  # vitesse du son en m/s

# Micros horizontaux : sur x, y=0
Coordonné_micro_h = 0.55 - 0.16 * np.arange(8)
y2 = np.zeros(8)  # y = 0 pour tous

# Micros verticaux : sur z (axe vertical), x=0, y=0
Coordonné_micro_v = 0.27 + 0.16 * np.arange(8)

# Optionnel : pour compatibilité 2D
y1 = Coordonné_micro_h  # ton code attend 'y1' pour former un array 'y'
# Matrice des coordonnées (8,2)
y = np.column_stack((y1, y2))

# Pour usage 3D croix
xm = y[:, 0]  # x des micros
ym = y[:, 1]  # y des micros

nb_points = 100
deltaXs = 1.12 / nb_points

indice_m = np.argmax(np.fft[:50000,0])

# On défini la fonction de Green
def G2(x, f, c):
    return(np.exp(-1j*(2*pi*f/c)*x)/(4*pi*x) )

xs = np.arange(nb_points) * deltaXs
ys = np.arange(nb_points) * deltaXs + 0.3
xs, yS = np.meshgrid(xs,ys)
M = np.fft[indice_m]
deltaX = 0.16
yl = np.arange(8)*deltaX
y2 = np.zeros([8,1])
y = np.column_stack((y1, y2))
xm = y[:,0]
ym = y[:,1]
deltax = xs.flatten() - xm[:,np.newaxis ]
deltay = ys.flatten() - ym[:,np.newaxis ]
N = np.sqrt(deltax**2 + deltay**2)
g = G2(N, 2000, c)
print(g.shape)

ng = np.sqrt(np.sum(np.abs(g)**2,axis = 0))
gn = g/ng

Z = (np.conj(gn).T) @ M
print(Z)
Z = np.reshape(Z, xs.shape)

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d' )
ax.plot_surface(xs, ys, np.abs(Z), cmap='viridis' )
# Ajout des titres aux axes
ax.set_xlabel('Axe X')
ax.set_ylabel('Axe Y')
ax.set_zlabel('vraissemblance')

plt.show()




X = np.load(r".\data\source2000.npy")
fft = np.fft.fft(X, axis=0)
a=1
nb_points = 100
deltaXs = 1.12/nb_points
deltaX = 0.16
y1 = np.arange(8)*deltax
y2 = np.zeros([8,1])
y = np.column_stack((y1, y2))
fe = 20000
c = 340
frequencies = np.fft.fftfreq(len(X[:,0]), 1/fe)
indice_m = np.argmax(fft[:50000,0])
f = frequencies[indice_m]

M = fft[indice_m]

f = 2000

Xsmv = np.zeros([nb_points,nb_points])
Y = [i*deltaXs for i in range(nb_points) ]
x = [i*deltaXs for i in range(nb_points) ]
x, Y = np.meshgrid(x, Y)
for i in range(nb_points-1):
    for j in range(nb_points-1):
        g = G2(y-[x[i,j], Y[i,j]], f, c)
        norme_g = np.linalg.norm(M@np.conj(g.T)/np.linalg.norm(g) )
        Xsmv[i,j] = norme_g

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d' )
ax.plot_surface(x, Y, Xsmv, cmap='viridis' )
ax.set_xlabel('Axe X')
ax.set_ylabel('Axe y')
ax.set_zlabel('vraissemblance' )
plt.show()




xcadrillage = np.arange(nb_points) * deltaXs - 0.6
ycadrillage = np.arange(nb_points) * deltaXs + 0.3

zcadrillage = np.arange(nb_points) * deltaXs + 1

xcadrillage, ycadrillage = np.meshgrid(xcadrillage, ycadrillage)
ycadrillage = np.arange(nb_points) * deltaXs + 0.3

zcadrillage, ycadrillage = np.meshgrid(zcadrillage, ycadrillage)

M = fft[indice_m]

deltax_croix_h = xcadrillage.flatten() - Coordonné_micro_h[:,np.newaxis ]
deltay_croix_h = ycadrillage.flatten() - y2[:,np.newaxis ]

deltaz_croix_v = zcadrillage.flatten() - Coordonné_micro_v[:,np.newaxis ]
deltay_croix_v = ycadrillage.flatten() - y2[:,np.newaxis ]

N_croix_h = np.sqrt(deltax_croix_h**2 + deltay_croix_h**2)
g_croix_h = G2(N_croix_h, 2000, c)
N_croix_v = np.sqrt(deltaz_croix_v**2 + deltay_croix_v**2)
g_croix_v = G2(N_croix_v, 2000, c)

ng_croix_h = np.sqrt(np.sum(np.abs(g_croix_h)**2,axis = 0))
gn_croix_h = g_croix_h/ng_croix_h

ng_croix_v = np.sqrt(np.sum(np.abs(g_croix_v)**2,axis = 0))
gn_croix_v = g_croix_v/ng_croix_v

Z_croix_h = (np.conj(gn_croix_h).T) @ M
Z_croix_h = np.reshape(Z_croix_h, xs.shape)
coord_h = np.argmax(Z_croix_h)

Z_croix_v = (np.conj(gn_croix_v).T) @ M
Z_croix_v = np.reshape(Z_croix_v, xs.shape)

coord_v = np.argmax(Z_croix_v)
coord_h = np.argmax(Z_croix_h)
coord_v = [coord_v//nb_points,coord_v%nb_points ]
coord_h = [coord_h//nb_points,coord_h%nb_points ]

print("les coordonnées de la source sont:", xcadrillage[coord_h[0], 0], ycadrillage[coord_v[1],0], zcadrillage[ coord_v[0], 0])
fig = plt.figure()
ax = fig.add_subplot(122, projection='3d' )
ax.plot_surface(xcadrillage, ycadrillage, np.abs(Z_croix_h), cmap='viridis' )
ax.set_xlabel('Axe X')
ax.set_ylabel('Axe Y')
ax.set_zlabel('vraissemblance')

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.plot_surface(xcadrillage, ycadrillage, np.abs(Z_croix_v), cmap='viridis' )
ax.set_xlabel('Axe Z')
ax.set_ylabel('Axe Y')
ax.set_zlabel('vraissemblance')
plt.show()

print("les coordonnées de la source sont:", xcadrillage[coord_h[0], 0], ycadrillage[coord_v[1], 0], zcadrillage[coord_v[0], 0])
