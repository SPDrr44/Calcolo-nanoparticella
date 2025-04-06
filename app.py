import streamlit as st
import numpy as np
from scipy.optimize import brentq

h = 6.62607015e-34
hbar = h / (2.0 * np.pi)
c = 3.0e8
e_charge = 1.602176634e-19
eps0 = 8.8541878128e-12
m0 = 9.10938356e-31
eV_to_J = e_charge

Eg_bulk_eV = 3.4
Eg_bulk_J  = Eg_bulk_eV * eV_to_J

me_eff = 0.24 * m0
mh_eff = 0.59 * m0
epsilon_r = 8.66

R_min = 0.1e-9
R_max = 10e-9

def Eg_confined_J(R):
    kinetic = (hbar**2 * (np.pi**2) / (2 * R**2)) * (1/me_eff + 1/mh_eff)
    coulomb = (1.8 * (e_charge**2)) / (4 * np.pi * eps0 * epsilon_r * R)
    return Eg_bulk_J + kinetic - coulomb

def solve_for_radius(Eg_target_J):
    def f(R): return Eg_confined_J(R) - Eg_target_J
    return brentq(f, R_min, R_max)

st.title("Calcolo Raggio Nanoparticella (ZnO)")
lambda_input = st.number_input("Lunghezza d'onda (nm)", value=350.0, min_value=0.0)

if st.button("Calcola"):
    lambda_m = lambda_input * 1e-9
    Eg_target = (h * c) / lambda_m

    try:
        R_m = solve_for_radius(Eg_target)
        R_nm = R_m * 1e9
        Eg_eV = Eg_confined_J(R_m) / eV_to_J

        st.write(f"**Raggio**: {R_nm:.3f} nm")
        st.write(f"**Gap calcolato**: {Eg_eV:.3f} eV")
    except Exception as e:
        st.error(f"Errore: {e}")
