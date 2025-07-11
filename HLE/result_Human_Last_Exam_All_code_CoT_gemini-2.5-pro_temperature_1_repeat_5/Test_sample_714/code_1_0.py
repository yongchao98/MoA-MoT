import numpy as np

# --- Constants ---
q = 1.602e-19  # Electron charge in C
a0 = 5.292e-11 # Bohr radius in m
h_bar = 1.0546e-34 # Reduced Planck constant in J*s
c = 2.998e8   # Speed of light in m/s
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m
lam = 589e-9      # Wavelength in m
tau_exp = 16.2e-9   # Experimental lifetime in s

# --- Problem Parameters ---
# For the 3p -> 3s transition in Sodium
n = 3
l_upper = 1 # 3p state

# We use an effective nuclear charge Z=1 for the valence electron,
# as it is screened by the 10 inner-shell electrons.
Z = 1.0

# --- Calculations ---

# 1. Angular frequency of the transition
omega = 2 * np.pi * c / lam

# 2. Calculate the squared radial integral |I_rad|^2
# For a hydrogenic transition between (n,l) and (n, l-1), the radial integral is:
# |<n,l-1|r|n,l>|^2 = ((3*n/2) * sqrt(n^2 - l^2) * a0/Z)^2
# For 3p -> 3s, n=3, l=1
I_rad_sq = ((3 * n / 2) * np.sqrt(n**2 - l_upper**2) * a0 / Z)**2

# 3. Calculate the averaged squared transition dipole moment |r_avg|^2
# |r_avg|^2 = (1/(2*l_upper+1)) * |I_rad|^2
g_upper = 2 * l_upper + 1
r_avg_sq = (1 / g_upper) * I_rad_sq

# 4. Calculate the spontaneous emission rate (Einstein A coefficient)
# A21 = (q^2 * omega^3) / (3 * pi * epsilon_0 * h_bar * c^3) * |r_avg|^2
A21_prefactor = (q**2 * omega**3) / (3 * np.pi * epsilon_0 * h_bar * c**3)
A21 = A21_prefactor * r_avg_sq

# 5. Calculate the theoretical lifetime
tau_theo = 1 / A21

# 6. Compare with the experimental lifetime
ratio = tau_theo / tau_exp

# --- Output the results ---
print("--- Calculation Steps ---")
print(f"1. Squared radial integral |I_rad|^2 = {I_rad_sq:.3e} m^2")
print(f"2. Averaged squared dipole moment |r_avg|^2 = {r_avg_sq:.3e} m^2")
print(f"3. Spontaneous emission rate A_21 = {A21:.3e} s^-1")
print("\n--- Final Result ---")
print(f"Theoretical lifetime (τ_theo): {tau_theo*1e9:.3f} ns")
print(f"Experimental lifetime (τ_exp): {tau_exp*1e9:.3f} ns")
print("\n--- Comparison ---")
print(f"The ratio of theoretical to experimental lifetime is:")
print(f"τ_theo / τ_exp = {tau_theo*1e9:.3f} ns / {tau_exp*1e9:.3f} ns = {ratio:.3f}")

print("\nThis means the calculated theoretical lifetime is approximately one-tenth as long as the experimentally measured lifetime.")
<<<C>>>