import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd

# experimental data, Si photodiode
T = np.log10(np.array([10, 30, 50, 100, 150, 200, 300])) # unit: K
df = pd.read_csv("data/Si.csv")

Capacitance_10K = df['10'].to_numpy() # unit: pF
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
eps = 11.7 # dielectric constant
me = 0.2*elecmass # electron effective mass
mh = 0.49*elecmass # hole effective mass
EG = 1.12*q # band gap in joules

def C_T_approx(T, a, b, c1, c2, c3):
    return a / np.sqrt((c1*T + c2*T*np.log(T) + c3)*(1 + b**(1/T)))

fig, ax = plt.subplots()

data_list = list([Capacitance_10K, Capacitance_30K, Capacitance_50K, Capacitance_100K, 
                  Capacitance_150K, Capacitance_200K, Capacitance_300K])

Capacitance_0V = np.array([elem[0] for elem in data_list]) # unit: pF
Capacitance_n1V = np.array([elem[1] for elem in data_list])
Capacitance_n2V = np.array([elem[2] for elem in data_list])
Capacitance_n3V = np.array([elem[3] for elem in data_list])
Capacitance_n4V = np.array([elem[4] for elem in data_list])

ax.scatter(T, Capacitance_0V*pF)
ax.scatter(T, Capacitance_n1V*pF)
ax.scatter(T, Capacitance_n2V*pF)
ax.scatter(T, Capacitance_n3V*pF)
ax.scatter(T, Capacitance_n4V*pF)

ax.legend(['0V', '-1V', '-2V', '-3V', '-4V'], loc='upper left')
ax.set_xlabel("log Temperature (K)")
ax.set_ylabel("Capacitance (F)")
ax.set_title("Measured capacitance of Si photodiode" "\n" 
             "at different temperatures from S22 measurement", loc='center')

plt.show()