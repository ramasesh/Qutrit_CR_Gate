import numpy as np
import qutip.three_level_atom as qutrit_utils
import qutip as qt

basis_states = qutrit_utils.three_level_basis()
eq_sup = 1/np.sqrt(3)*(basis_states[0] + basis_states[1] + basis_states[2])
sig00, sig11, sig22, sig01, sig12 = qutrit_utils.three_level_ops()
sig02 = qt.Qobj(np.array([[0,0,1],[0,0,0],[0,0,0]])) 
sig02x = sig02 + sig02.dag()
sig02y = +1j*sig02 - 1j*sig02.dag()
sig01x = sig01 + sig01.dag()
sig01y = +1j*sig01 - 1j*sig01.dag()
sig12x = sig12 + sig12.dag()
sig12y = +1j*sig12 - 1j*sig12.dag()
