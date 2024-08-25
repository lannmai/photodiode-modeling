import numpy as np
import matplotlib.pyplot as plt
import scipy 

# experimental data, Si photodiode

V = np.array([0, -1, -2, -3, -4]) # unit: V
Capacitance_6K = np.array([1.36, 1.29, 1.26, 1.23, 1.2]) # unit: pF
Capacitance_10K = np.array([1.37, 1.3, 1.26, 1.24, 1.21])
Capacitance_30K = np.array([1.51, 1.41, 1.37, 1.33, 1.29])
Capacitance_77K = np.array([2.68, 1.92, 1.82, 1.76, 1.71])
Capacitance_300K = np.array([11.8, 2.73, 2.36, 2.17, 2.06])

# physical constants
elecmass = scipy.constants.electron_mass
q = scipy.constants.elementary_charge
k = scipy.constants.Boltzmann
h = scipy.constants.Planck
pF = 1e-12

# material properties
eps = 11.7 # dielectric constant
me = 0.2*elecmass # electron effective mass
mh = 0.49*elecmass # hole effective mass
EG = 1.12*q # band gap in joules

def C_total(V, a, V_0, C_p):
    return a/np.sqrt(V_0-V)+C_p

data_list = list([Capacitance_6K, Capacitance_10K, Capacitance_30K, Capacitance_77K, Capacitance_300K])
fit_params = list()
fit_params_pcov = list()

for data in data_list:
    fit_params.append(scipy.optimize.curve_fit(C_total, V, data*1e-12)[0]) # fitted params
    fit_params_pcov.append(scipy.optimize.curve_fit(C_total, V, data*pF)[1]) # fitted params' covariances

fig, ax = plt.subplots()

V_fine = np.linspace(min(V), max(V), num=400) # finer voltage step for plotting model

ax.scatter(V, Capacitance_6K*pF, c='red')
ax.scatter(V, Capacitance_10K*pF, c='green')
ax.scatter(V, Capacitance_30K*pF, c='blue')
ax.scatter(V, Capacitance_77K*pF, c='cyan')
ax.scatter(V, Capacitance_300K*pF, c='magenta')

ax.plot(V_fine, C_total(V_fine, fit_params[0][0], fit_params[0][1], fit_params[0][2]), 'r--')
ax.plot(V_fine, C_total(V_fine, fit_params[1][0], fit_params[1][1], fit_params[1][2]), 'g--')
ax.plot(V_fine, C_total(V_fine, fit_params[2][0], fit_params[2][1], fit_params[2][2]), 'b--')
ax.plot(V_fine, C_total(V_fine, fit_params[3][0], fit_params[3][1], fit_params[3][2]), 'c--')
ax.plot(V_fine, C_total(V_fine, fit_params[4][0], fit_params[4][1], fit_params[4][2]), 'm--')

ax.legend(['6K', '10K', '30K', '77K', '300K'], loc='upper left')
ax.set_xlabel("Applied reverse bias (V)")
ax.set_ylabel("Capacitance (F)")
ax.set_title("Measured (scatter points) and modeled (dashed lines) capacitance of" "\n" 
             "Si photodiode at different temperatures from S22 measurement", loc='center')
plt.show()
