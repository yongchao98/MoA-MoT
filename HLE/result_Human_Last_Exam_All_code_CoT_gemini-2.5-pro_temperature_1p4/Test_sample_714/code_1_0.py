import numpy as np

# Step 1: Define constants and parameters
q = 1.602e-19  # Electron charge in C
a0 = 5.292e-11 # Bohr radius in m
hbar = 1.054e-34 # Reduced Planck constant in J.s
c = 3.00e8     # Speed of light in m/s
eps0 = 8.854e-12 # Permittivity of free space in F/m
Z = 1.0        # Effective nuclear charge for Sodium valence electron

# Problem-specific parameters
lambda_nm = 589  # Wavelength in nm
lambda_m = lambda_nm * 1e-9  # Wavelength in m
tau_exp_ns = 16.2 # Experimental lifetime in ns
tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in s

# Quantum numbers for the transition 3p -> 3s
n1, l1 = 3, 0  # Lower state (3s)
n2, l2 = 3, 1  # Upper state (3p)

# Step 2: Calculate the radial integral squared
# The analytical result for the radial integral for n,l -> n,l-1 is (3/2Z) * n * sqrt(n^2 - l^2) * a0
# For 3p -> 3s, n=3, l=1
radial_integral = (3/2.0 * n2 * np.sqrt(n2**2 - l2**2)) * (a0 / Z)
radial_integral_sq = radial_integral**2

# Step 3: Calculate the Einstein A coefficient
# Angular frequency
omega = 2 * np.pi * c / lambda_m

# Degeneracy of the upper state (ignoring spin) and max orbital quantum number
g2 = 2 * l2 + 1
l_max = max(l1, l2)

# Einstein A coefficient A_21
# A_21 = (e^2 * omega^3) / (3 * pi * eps0 * hbar * c^3) * (l_max/g_2) * I_r^2
A_21 = (q**2 * omega**3) / (3 * np.pi * eps0 * hbar * c**3) * (l_max / g2) * radial_integral_sq

# Step 4: Calculate the theoretical lifetime
tau_theo_s = 1 / A_21
tau_theo_ns = tau_theo_s * 1e9

# Step 5: Compare theoretical and experimental lifetimes
ratio = tau_theo_s / tau_exp_s

# Print the results
print("Calculation Steps and Results:")
print(f"Transition: {n2}p -> {n1}s")
print(f"Effective Nuclear Charge (Z_eff) used: {Z}")
print(f"Radial Integral (I_r): {radial_integral:.4e} m")
print(f"Squared Radial Integral (I_r^2): {radial_integral_sq:.4e} m^2")
print(f"Angular Frequency (omega): {omega:.4e} rad/s")
print(f"Einstein A Coefficient (A_21): {A_21:.4e} s^-1")
print("-" * 30)
print(f"Theoretical Lifetime (calculated): {tau_theo_ns:.2f} ns")
print(f"Experimental Lifetime (given): {tau_exp_ns:.2f} ns")
print("-" * 30)
print("Comparison:")
print(f"The ratio (Theoretical / Experimental) is {tau_theo_ns:.2f} ns / {tau_exp_ns:.2f} ns = {ratio:.3f}")
print("\nThe calculated theoretical lifetime is approximately 0.1 times the experimental lifetime.")
<<<C>>>