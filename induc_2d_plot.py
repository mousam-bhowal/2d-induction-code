import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Change the filenames here, as given in induc_2d.f90

data = np.loadtxt("induc_2d_dat_1")
x = np.loadtxt("induc_2d_Xvals_1")
y = np.loadtxt("induc_2d_Yvals_1")
eta_t = np.loadtxt("induc_2d_eta_t_1")

filename = "induc_2d_plot_1.gif" # this is the output file, change the name appropriately for a new plot.

# Change the title of the plot, which contains the values of velocity field and diffusivity.
title = r'Plot of $T(t, x, y), \quad \mathbf{v} = 5\hat{x} + \hat{y},\quad \eta = 0.3$'

eta = eta_t[0]
t = eta_t[1:]


X, Y = np.meshgrid(x, y)
Nx, Ny, steps = len(x), len(y), len(t)


T = np.zeros((steps, Nx, Ny))

for n in range(steps):
	for i in range(Nx):
		T[n][i] = data[Nx*n + i]

fig = plt.figure(figsize = (6, 5))
ax = fig.add_subplot(111, projection = '3d')

ax.view_init(elev = 15, azim = -45)
ax.set_zlim(0, 10.1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$T$')


ax.set_title(title)


lvl = 50
contour = ax.contour3D(X, Y, T[0].T, levels = lvl)
text = ax.text(min(x), min(y), 9.5, r'$t = $%.3f' % 0)


def update(frame):
	global contour
	global text
	contour.remove()
	text.remove()
	
	contour = ax.contour3D(X, Y, T[frame].T, levels = lvl)
	text = ax.text(min(x), min(y), 9.5, r'$t = $%.2f' % t[frame])
#	return contour.collections

fig.tight_layout()
ani = FuncAnimation(fig, update, frames = range(steps), interval = 500, blit = False)
ani.save(filename, writer = 'pillow', fps = 2, dpi = 120)

plt.tight_layout()
plt.show()
