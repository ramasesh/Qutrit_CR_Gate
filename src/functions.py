# Contains the functions necessary for the time-dependent Hamiltonian simulations

import numpy as np

def cos_pulse_functions(t_pulse=0.3, epsilon_x=10):

    def eps_x(t):
        return -1 * (np.cos(2*np.pi*t/t_pulse) - 1) * (np.heaviside(t,1) - np.heaviside(t - t_pulse,1)) * epsilon_x
    def deps_x(t):
        return 2*np.pi/t_pulse * (np.sin(2*np.pi*t/t_pulse)) * (np.heaviside(t,1) - np.heaviside(t - t_pulse,1)) * epsilon_x
    def eps_y(t):
        return 0

    return [eps_x, deps_x, eps_y]
