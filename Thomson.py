import numpy as np
from random import random as rnd
from scipy.spatial import SphericalVoronoi
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

N0 = 188

def init(N):
    points = np.zeros((N,2))
    for i in range(N):
        points[i, 0] = -np.pi/2.0 + np.pi*rnd()
        points[i, 1] = 2*np.pi*rnd()
    return points

def phi_pair(angs1, angs2):
    dx = np.sin(angs1[0])*np.cos(angs1[1]) - np.sin(angs2[0])*np.cos(angs2[1])
    dy = np.sin(angs1[0])*np.sin(angs1[1]) - np.sin(angs2[0])*np.sin(angs2[1])
    dz = np.cos(angs1[0]) - np.cos(angs2[0])
    return 1.0/np.sqrt(dx**2 + dy**2 + dz**2)

def phi1(points):
    phi = 0.0
    for i in range(len(points)):
        for j in range(len(points)):
            if i < j:
                phi += phi_pair(points[i], points[j])
    return phi
def fill_pairwise(pairwise, points):
    for i in range(len(points)):
        for j in range(len(points)):
            if i < j:
                pairwise[i][j] = phi_pair(points[i], points[j])
    return pairwise
def change_i_pairs(i, phi_pairwise, points):
    for j in range(len(points)):
        if i < j:
            phi_pairwise[i][j] = phi_pair(points[i], points[j])
        if i > j:
            phi_pairwise[j][i] = phi_pair(points[i], points[j])
    return phi_pairwise

def derive(points, i, ang):
    der = 0
    pi0 = points[i][0]
    pi1 = points[i][1]
    if ang == 0: #theta
        for j in xrange(N0):
            if i != j:
                pj0 = points[j][0]
                pj1 = points[j][1]
                dx = np.sin(pi0)*np.cos(pi1) - np.sin(pj0)*np.cos(pj1)
                dy = np.sin(pi0)*np.sin(pi1) - np.sin(pj0)*np.sin(pj1)
                dz = np.cos(pi0) - np.cos(pj0)
                pair = -1.0/np.sqrt(dx**2 + dy**2 + dz**2)**3
                pair *= dx*np.cos(pi0)*np.cos(pi1)+dy*np.cos(pi0)*np.sin(pi1)-dz*np.sin(pi0)
                der += pair
    elif ang == 1: #phi
        for j in xrange(N0):
            if i != j:
                pj0 = points[j][0]
                pj1 = points[j][1]
                dx = np.sin(pi0)*np.cos(pi1) - np.sin(pj0)*np.cos(pj1)
                dy = np.sin(pi0)*np.sin(pi1) - np.sin(pj0)*np.sin(pj1)
                dz = np.cos(pi0) - np.cos(pj0)
                pair = -1.0/np.sqrt(dx**2 + dy**2 + dz**2)**3
                pair *= -dx*np.sin(pi0)*np.sin(pi1)+dy*np.sin(pi0)*np.cos(pi1)
                der += pair
    return der

def myoptimizer(N, maxiter = 200, points = [], alpha = 0.003, maxsubiter = 1):
    #Gets N - number of points and computes minimal energy with coordinate descent algorithm
    #1-D optimizer is gradient descent algorithm with parameter alpha = 0.5,
    #which is decreasing over iterations
    #maxsubiter = 100 defines the number of iterations for 1-D optimizaer
    if points == []:
        points = init(N)
    phi_old = phi1(points)
    phi_new = phi_old
    #phi_pairwise = fill_pairwise(np.zeros((N,N)), points)
    energy = []
    energy.append(phi_new)
    for it in xrange(maxiter):
        for i in xrange(len(points)):
            for subit in xrange(maxsubiter):
                der0 = derive(points, i, 0)
                der1 = derive(points, i, 1)
                points[i][0] -= alpha*der0#/(1+it)
                points[i][1] -= alpha*der1#/(1+it)
        #phi_pairwise = change_i_pairs(i, phi_pairwise, points)
        #print phi_pairwise
        #phi_new = phi_pairwise.sum()
        phi_new = phi1(points)
        energy.append(phi_new)
        if (np.abs(phi_new - phi_old) < 1e-80): break
        else: phi_old = phi_new
    return phi_new, points, energy#, phi_pairwise


#alphas = np.linspace(0.008, 0.1, 300)
#phis = []
polarpoints = np.zeros((N0,2))
for i in range(N0):
    polarpoints[i, 0] = 0.2*rnd()
    polarpoints[i, 1] = 2*np.pi*rnd()
import time
start = time.time()

#for alpha in alphas:
phi, points, energy = myoptimizer(N = N0, points = polarpoints)
    #phis.append(phi)

end = time.time()
print 't:', end-start

print phi
print points

'''plt.plot(alphas, phis)
plt.title('constant step (different alphas)')
plt.xlabel('alpha')
plt.ylabel('min energy')
plt.grid(True)
plt.show()
'''
import pickle
with open('points constant step alpha: 0003, '+str(N0)+'.pk', 'wb') as f:
    pickle.dump(points, f)

fig, axs = plt.subplots(2, 1)
#axs[0].title('fsd')
axs[0].plot(xrange(len(energy)), energy)

axs[0].set_ylabel('energy')
axs[1].plot(xrange(3, len(energy)), energy[3:])
axs[1].set_xlabel('iterations')
axs[1].set_ylabel('energy')
axs[0].grid(True)
axs[1].grid(True)
axs[0].title.set_text('constant step alpha: 0.3')
axs[1].title.set_text('zoom in')
plt.savefig('N:' + str(N0) + ' Energy-iterations graphs (constant step alpha: 0.003).jpg', dpi = 200)
plt.show()

pointsXYZ = np.zeros((N0, 3))
for i in range(len(points)):
    X = np.sin(points[i][0])*np.cos(points[i][1])
    Y = np.sin(points[i][0])*np.sin(points[i][1])
    Z = np.cos(points[i][0])
    pointsXYZ[i][0] = X
    pointsXYZ[i][1] = Y
    pointsXYZ[i][2] = Z

sv = SphericalVoronoi(pointsXYZ)

sv.sort_vertices_of_regions()
# generate plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# plot the unit sphere for reference (optional)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='y', alpha=0.1)
# plot generator points
ax.scatter(pointsXYZ[:, 0], pointsXYZ[:, 1], pointsXYZ[:, 2], c='b', s = 2)
# plot Voronoi vertices
#ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')
# indicate Voronoi regions (as Euclidean polygons)
for region in sv.regions:
    random_color = colors.rgb2hex(np.random.rand(3))
    polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
    if len(sv.vertices[region]) == 3:
        polygon.set_color('pink')
    elif len(sv.vertices[region]) == 4:
        polygon.set_color('yellow')
    elif len(sv.vertices[region]) == 5:
        polygon.set_color('red')
    elif len(sv.vertices[region]) == 6:
        polygon.set_color('green')
    elif len(sv.vertices[region]) == 7:
        polygon.set_color('black')
    else:
        polygon.set_color(random_color)
    polygon.set_edgecolor('black')
    ax.add_collection3d(polygon)
plt.show()
plt.savefig('N0 ' + str(N0)+ ' Voronoi (constant step alpha 0003)', dpi = 200)
