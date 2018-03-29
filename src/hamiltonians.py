import sys
import os
import numpy as np
import qutip as qt
from src import qutrit 

def H0(delta_1, delta_2, lambda_t):
    """ Gives the static part of Hamiltonian in the rotating frame for a 
    transmon driven by H_d = epsilon_x ( sigma_x_01 + lambda sigma_x_12 ) """
    
    return delta_1*qutrit.sig11 + delta_2*qutrit.sig22

def H1(lambda_t):
    """ Gives the dynamic part of the Hamiltonian in the rotating 
    frame for a transmon driven by H_d = epsilon_x (sigma_x_01 + lambda sigma_x_12 ) """
    
    return (qutrit.sig01x + lambda_t * qutrit.sig12x)/2. 

def H2(lambda_t):
    """ Gives the dynamic part of the Hamiltonian in the rotating 
    frame for a transmon driven by H_d = epsilon_x (sigma_y_01 + lambda sigma_y_12 ) """
    
    return (qutrit.sig01y + lambda_t * qutrit.sig12y)/2. 

def HR(eps_x_func, deps_x_func, eps_y_func, alpha, delta_1, delta_2, lambda_t):
    """ returns an array of Hamiltonians, with their time-dependent coefficients,
    which represent the initial (untransformed) qutrit Hamiltonian """

    def coeff_1(t, args):
        return eps_x_func(t)
    def coeff_2(t, args):
        return eps_y_func(t)

    return [H0(delta_1, delta_2, lambda_t),[H1(lambda_t), coeff_1], [H2(lambda_t), coeff_2]]

def HV(eps_x_func, deps_x_func, eps_y_func, alpha, delta_1, delta_2, lambda_t):
    """ returns an array of Hamiltonians, with their time-dependent coefficients,
    which represent the DRAG frame transformed qutrit Hamiltonian """
    
    HV_1 = qutrit.sig01x
    HV_2 = qutrit.sig02x
    HV_3 = qutrit.sig22 
    HV_4 = qutrit.sig11 
    HV_5 = (qutrit.sig01y + lambda_t*qutrit.sig12y) 

    def coeff_1(t, args):
        return eps_x_func(t)/2.
    def coeff_2(t, args):
        return lambda_t*np.power(eps_x_func(t),2)/(8*alpha)
    def coeff_3(t, args):
        return (delta_2 + (np.power(lambda_t,2) + 2)/(4*alpha)*np.power(eps_x_func(t),2))
    def coeff_4(t, args):
        return (delta_1 - (np.power(lambda_t,2) - 4)/(4*alpha)*np.power(eps_x_func(t),2))
    def coeff_5(t, args):
        return (eps_y_func(t)/2. + deps_x_func(t)/(2*alpha))

    return [[HV_1, coeff_1], [HV_2, coeff_2], [HV_3, coeff_3], [HV_4, coeff_4] , [HV_5, coeff_5]]

def H_transmon(omega_01, alpha, n_levels):
    """ returns the bare transmon Hamiltonian """ 

    a = qt.destroy(n_levels)
    H = omega_01 * a.dag() * a + alpha * a.dag() * a.dag() * a * a

    return [H, a]

def H_two_transmons(omega_1, alpha_1, omega_2, alpha_2, g, n_levels):
    """ returns the Hamiltonian for two coupled transmons """

    H1, a1 = H_transmon(omega_1, alpha_1, n_levels)
    H2, a2 = H_transmon(omega_2, alpha_2, n_levels)
   
    H1 = qt.tensor(qt.identity(n_levels), H1)
    H2 = qt.tensor(H2, qt.identity(n_levels))
    a1 = qt.tensor(qt.identity(n_levels), a1)
    a2 = qt.tensor(a2, qt.identity(n_levels))

    H_total = H1 + H2 + g*(a1.dag() * a2 + a2.dag() * a1)

    return [H_total, a1, a2]

def transmon_drive(drive_function, dest_operator):
    """ given a destruction operator and a time-dependent drive function, constructs the time-dependent 
    Hamiltonian representing the drive """
    
    return [dest_operator + dest_operator.dag(), drive_function]

def H_two_transmons_with_drive(omega_1, alpha_1, omega_2, alpha_2, g, n_levels, drive_function):
    """ Constructs the full Hamiltonian for two transmons, with a drive on the first one """ 
    [H_static, a1, a2] = H_two_transmons(omega_1, alpha_1, omega_2, alpha_2, g, n_levels) 
    H_drive = transmon_drive(drive_function, a1)

    return [H_static, H_drive]
