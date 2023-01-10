# First set up the figure, the axis, and the plot element we want to animate
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

fig = plt.figure()
ax = plt.axes(xlim=(0, 10), ylim=(0, 10))
# line, = ax.plot([], [], lw=2)
a = np.random.random((5, 5))
b = np.random.random((10, 10))

im = [plt.imshow(a, interpolation='none'), plt.imshow(b, interpolation='none')]


# initialization function: plot the background of each frame
# initialization function: plot the background of each frame
def init():
    im[0].set_data(np.random.random((5, 5)))

    return [im]


# animation function.  This is called sequentially
def animate(i):
    print(i)
    a = im[0].get_array()
    a = a * np.exp(-0.001 * i)  # exponential decay of the values
    im[0].set_array(a)

    im[1].set_data(np.random.random((10, 2)))
    return [im]


init()
ani = FuncAnimation(fig, animate, save_count=100)

plt.show()
