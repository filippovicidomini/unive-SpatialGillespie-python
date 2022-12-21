# i want to simulae the reaction of belousov zhabotinsky reaction in python.
# and show it in a 2D plot.

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

# define the size of the grid
size = 100


# define the initial state of the grid
def initial_state():
    state = np.zeros((size, size))
    state[50:60, 50:60] = 2
    state[30:35, 30:35] = 2
    return state


# define the update function
def update(state):
    new_state = state.copy()
    for i in range(size):
        for j in range(size):
            total = (state[i, (j - 1) % size] + state[i, (j + 1) % size] +
                     state[(i - 1) % size, j] + state[(i + 1) % size, j])
            if total == 0:
                new_state[i, j] = 0
            elif total == 1:
                new_state[i, j] = 1
            elif total == 2:
                new_state[i, j] = 2
            elif total == 3:
                new_state[i, j] = 0
            elif total == 4:
                new_state[i, j] = 1
    return new_state


# define the animation function
def animate(i):
    global state
    state = update(state)
    im.set_data(state)
    return im,


# define the main function
def main():
    global state, im
    state = initial_state()
    fig, ax = plt.subplots()
    im = ax.imshow(state, interpolation='nearest', cmap='jet')
    ani = animation.FuncAnimation(fig, animate, interval=1, save_count=100)
    ani.save('bzr.gif', writer='imagemagick', fps=10)
    plt.show()


if __name__ == '__main__':
    main()
