import numpy as np
from random import random as rnd
from scipy.spatial import SphericalVoronoi
from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d

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

def derive(points, i, ang):
    der = 0
    if ang == 0: #theta
        for j in range(len(points)):
            if i != j:
                dx = np.sin(points[i][0])*np.cos(points[i][1]) - np.sin(points[j][0])*np.cos(points[j][1])
                dy = np.sin(points[i][0])*np.sin(points[i][1]) - np.sin(points[j][0])*np.sin(points[j][1])
                dz = np.cos(points[i][0]) - np.cos(points[j][0])
                pair = -1.0/np.sqrt(dx**2 + dy**2 + dz**2)**3
                pair *= dx*np.cos(points[i][0])*np.cos(points[i][1])+dy*np.cos(points[i][0])*np.sin(points[i][1])-dz*np.sin(points[i][0])
                der += pair
    elif ang == 1: #phi
        for j in range(len(points)):
            if i != j:
                dx = np.sin(points[i][0])*np.cos(points[i][1]) - np.sin(points[j][0])*np.cos(points[j][1])
                dy = np.sin(points[i][0])*np.sin(points[i][1]) - np.sin(points[j][0])*np.sin(points[j][1])
                dz = np.cos(points[i][0]) - np.cos(points[j][0])
                pair = -1.0/np.sqrt(dx**2 + dy**2 + dz**2)**3
                pair *= -dx*np.sin(points[i][0])*np.sin(points[i][1])+dy*np.sin(points[i][0])*np.cos(points[i][1])
                der += pair
    return der

def myoptimizer(N, maxiter = 200, points = [], alpha = 0.3, maxsubiter = 1):
    #Gets N - number of points and computes minimal energy with coordinate descent algorithm
    #1-D optimizer is gradient descent algorithm with parameter alpha = 0.5,
    #which is decreasing over iterations
    #maxsubiter = 100 defines the number of iterations for 1-D optimizaer

    if points == []:
        points = init(N)
    phi_old = phi1(points)
    phi_new = phi_old
    error = []
    error.append(phi_new)
    for it in range(maxiter):
        for i in range(len(points)):
            for subit in range(maxsubiter):
                der0 = derive(points, i, 0)
                der1 = derive(points, i, 1)
                points[i][0] -= alpha*der0/(1+it)
                points[i][1] -= alpha*der1/(1+it)
        phi_new = phi1(points)
        error.append(phi_new)
        if False*((phi_new - phi_old) < 1e-18): break
        else: phi_old = phi_new
    return phi_new, points, error

N0 = 327

polarpoints = np.zeros((N0,2))
for i in range(N0):
    polarpoints[i, 0] = 0.2*rnd()
    polarpoints[i, 1] = 2*np.pi*rnd()
import time
start = time.time()
phi, points, error = myoptimizer(N = N0, points = polarpoints)
end = time.time()
print 't:', end-start
#phi = 1
#points = polarpoints
print phi
print points

fig, axs = plt.subplots(2, 1)
#axs[0].title('fsd')
axs[0].plot(xrange(len(error)), error)

axs[0].set_ylabel('energy')
axs[1].plot(xrange(3, len(error)), error[3:])
axs[1].set_xlabel('iterations')
axs[1].set_ylabel('energy')
axs[0].grid(True)
axs[1].grid(True)
axs[0].title.set_text('Energy-iterations graphs (constant step)')
axs[1].title.set_text('zoom in (starting from 1 iteration)')
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
ax.scatter(pointsXYZ[:, 0], pointsXYZ[:, 1], pointsXYZ[:, 2], c='b')
# plot Voronoi vertices
ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')
# indicate Voronoi regions (as Euclidean polygons)
for region in sv.regions:
    random_color = colors.rgb2hex(np.random.rand(3))
    polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
    if len(sv.vertices[region]) == 3:
        polygon.set_color('yellow')
    elif len(sv.vertices[region]) == 4:
        polygon.set_color('green')
    elif len(sv.vertices[region]) == 5:
        polygon.set_color('blue')
    elif len(sv.vertices[region]) == 6:
        polygon.set_color('red')
    elif len(sv.vertices[region]) == 7:
        polygon.set_color('black')
    else:
        polygon.set_color(random_color)
    polygon.set_edgecolor('black')
    ax.add_collection3d(polygon)
plt.show()
