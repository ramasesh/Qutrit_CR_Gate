import sys
import os
import numpy as np
import qutip as qt
from src import qutrit 
from src import hamiltonians

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

    HR_t = hamiltonians.HR(eps_x_fun, deps_x_fun, eps_y_fun, alpha, delta_1, delta_2, lambda_t)
   
    HV_t = hamiltonians.HV(eps_x_fun, deps_x_fun, eps_y_fun, alpha, delta_1, delta_2, lambda_t)

    output1 = qt.mesolve(HR_t,psi0, t, [], [qutrit.sig00, qutrit.sig11, qutrit.sig22, qutrit.sig01x, qutrit.sig12x], options=options)

    output2 = qt.mesolve(HV_t, psi0, t, [],  [qutrit.sig00, qutrit.sig11, qutrit.sig22, qutrit.sig01x, qutrit.sig12x, qutrit.sig02x, qutrit.sig02y], options=options)

    return output1, output2
