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

def solve_dynamics(psi0=qutrit.eq_sup, t = np.linspace(-0.1,10.3,400), **kwargs): 
    omega = kwargs.pop('omega', 5000) #MHz
    alpha = kwargs.pop('alpha', -300) #MHz
    drive_omega = kwargs.pop('drive_omega', 5000) #MHz
    lambda_t = kwargs.pop('lambda_t', np.sqrt(2)) 

    eps_x_fun = kwargs['eps_x_fun']
    deps_x_fun = kwargs['deps_x_fun']
    eps_y_fun = kwargs['eps_y_fun']

    options = qt.Options()
    options.store_states = True

    delta_1 = omega - drive_omega
    delta_2 = alpha + 2 * delta_1

    HR_t = HR(eps_x_fun, deps_x_fun, eps_y_fun, alpha, delta_1, delta_2, lambda_t)
   
    HV_t = HV(eps_x_fun, deps_x_fun, eps_y_fun, alpha, delta_1, delta_2, lambda_t)

    output1 = qt.mesolve(HR_t,psi0, t, [], [qutrit.sig00, qutrit.sig11, qutrit.sig22, qutrit.sig01x, qutrit.sig12x], options=options)

    output2 = qt.mesolve(HV_t, psi0, t, [],  [qutrit.sig00, qutrit.sig11, qutrit.sig22, qutrit.sig01x, qutrit.sig12x, qutrit.sig02x, qutrit.sig02y], options=options)

    return output1, output2
