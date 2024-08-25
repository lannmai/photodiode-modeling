import numpy as np
import matplotlib.pyplot as plt
import scipy

# experimental data, InGaAs photodiode
V = np.array([0, -1, -2, -3, -4]) # unit: V
Capacitance_6K = np.array([5.84433, 5.72469, 5.65429, 5.59801, 5.55675]) # unit: pF
Capacitance_10K = np.array([5.7163, 5.6138, 5.56758, 5.52729, 5.4953])
Capacitance_30K = np.array([10.507, 9.3827, 8.43078, 7.74387, 7.2123])
Capacitance_50K = np.array([14.4979, 10.6152, 8.8966, 8.0028, 7.4594])
Capacitance_100K = np.array([19.1006, 12.0348, 9.72198, 8.61538, 7.9325])
Capacitance_300K = np.array([21.4766, 12.7451, 10.0689, 8.84538, 8.12815])

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

data_list = list([Capacitance_6K, Capacitance_10K, Capacitance_30K, Capacitance_50K, Capacitance_100K, Capacitance_300K])
fit_params = list()
fit_params_pcov = list()

for data in data_list:
    fit_params.append(scipy.optimize.curve_fit(C_total, V, data*1e-12)[0]) # fitted params
    fit_params_pcov.append(scipy.optimize.curve_fit(C_total, V, data*pF)[1]) # fitted params' covariances

fig, ax = plt.subplots()

V_fine = np.linspace(min(V), max(V), num=400) # finer voltage step for plotting model

C_6K = ax.scatter(V, Capacitance_6K*pF, c='red')
C_10K = ax.scatter(V, Capacitance_10K*pF, c='green')
C_30K = ax.scatter(V, Capacitance_30K*pF, c='blue')
C_50K = ax.scatter(V, Capacitance_50K*pF, c='cyan')
C_100K = ax.scatter(V, Capacitance_100K*pF, c='black')
C_300K = ax.scatter(V, Capacitance_300K*pF, c='magenta')

C_6K_model = ax.plot(V_fine, C_total(V_fine, fit_params[0][0], fit_params[0][1], fit_params[0][2]), 'r--')
C_10K_model = ax.plot(V_fine, C_total(V_fine, fit_params[1][0], fit_params[1][1], fit_params[1][2]), 'g--')
C_30K_model = ax.plot(V_fine, C_total(V_fine, fit_params[2][0], fit_params[2][1], fit_params[2][2]), 'b--')
C_50K_model = ax.plot(V_fine, C_total(V_fine, fit_params[3][0], fit_params[3][1], fit_params[3][2]), 'c--')
C_100K_model = ax.plot(V_fine, C_total(V_fine, fit_params[4][0], fit_params[4][1], fit_params[4][2]), 'b--')
C_300K_model = ax.plot(V_fine, C_total(V_fine, fit_params[5][0], fit_params[5][1], fit_params[5][2]), 'm--')

ax.legend(['6K', '10K', '30K', '50K', '100K', '300K'], loc='upper left')
ax.set_xlabel("Applied reverse bias (V)")
ax.set_ylabel("Capacitance (F)")
ax.set_title("Measured (scatter points) and modeled (dashed lines) capacitance of" "\n" 
             "InGaAs photodiode at different temperatures from S22 measurement", loc='center')
plt.show()