import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd

# experimental data, InGaAs photodiode
V = np.array([0, -1, -2, -3, -4]) # unit: V
df = pd.read_csv("data/InGaAs.csv")

Capacitance_6K = df['6'].to_numpy() # unit: pF
Capacitance_10K = df['10'].to_numpy()
Capacitance_30K = df['30'].to_numpy()
Capacitance_50K = df['50'].to_numpy()
Capacitance_100K = df['100'].to_numpy()
Capacitance_150K = df['150'].to_numpy()
Capacitance_200K = df['200'].to_numpy()
Capacitance_300K = df['300'].to_numpy()

# physical constants
elecmass = scipy.constants.electron_mass
q = scipy.constants.elementary_charge
k = scipy.constants.Boltzmann
h = scipy.constants.Planck
pF = 1e-12 # picofarad -> farad

# material properties
eps = 13.9 # dielectric constant
me = 0.043*elecmass # electron effective mass
mh = 0.46*elecmass # hole effective mass
EG = 0.74*q # band gap

def C_total(V, a, V_0, C_p):
    return a/np.sqrt(V_0-V)+C_p

data_list = list([Capacitance_6K, Capacitance_10K, Capacitance_30K, Capacitance_50K, 
                  Capacitance_100K, Capacitance_150K, Capacitance_200K, Capacitance_300K])
fit_params = list()
fit_params_pcov = list()

for data in data_list:
    fit_params.append(scipy.optimize.curve_fit(C_total, V, data*1e-12)[0]) # fitted params
    fit_params_pcov.append(scipy.optimize.curve_fit(C_total, V, data*pF)[1]) # fitted params' covariances

fig, ax = plt.subplots()

V_fine = np.linspace(min(V), max(V), num=400) # finer voltage step for plotting model

ax.scatter(V, Capacitance_6K*pF)
ax.scatter(V, Capacitance_10K*pF)
ax.scatter(V, Capacitance_30K*pF)
ax.scatter(V, Capacitance_50K*pF)
ax.scatter(V, Capacitance_100K*pF)
ax.scatter(V, Capacitance_150K*pF)
ax.scatter(V, Capacitance_200K*pF)
ax.scatter(V, Capacitance_300K*pF)

ax.plot(V_fine, C_total(V_fine, fit_params[0][0], fit_params[0][1], fit_params[0][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[1][0], fit_params[1][1], fit_params[1][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[2][0], fit_params[2][1], fit_params[2][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[3][0], fit_params[3][1], fit_params[3][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[4][0], fit_params[4][1], fit_params[4][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[5][0], fit_params[5][1], fit_params[5][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[6][0], fit_params[6][1], fit_params[6][2]))
ax.plot(V_fine, C_total(V_fine, fit_params[7][0], fit_params[7][1], fit_params[7][2]))

ax.legend(['6K', '10K', '30K', '50K', '100K', '150K', '200K', '300K'], loc='upper left')
ax.set_xlabel("Applied reverse bias (V)")
ax.set_ylabel("Capacitance (F)")
ax.set_title("Measured (scatter points) and modeled (dashed lines) capacitance of" "\n" 
             "InGaAs photodiode at different reverse bias from S22 measurement", loc='center')

plt.show()