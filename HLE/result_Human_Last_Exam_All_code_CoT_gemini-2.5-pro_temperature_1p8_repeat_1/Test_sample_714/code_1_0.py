import numpy as np

# Physical constants in SI units
h_bar = 1.054571817e-34  # J*s
c = 2.99792458e8         # m/s
epsilon_0 = 8.854187817e-12 # F/m
e = 1.602176634e-19      # C
a_0 = 5.29177210903e-11  # m (Bohr radius)

# Problem parameters
lambda_wav = 589e-9      # m
tau_exp = 16.2e-9        # s
Z = 11

# --- Step 1: Calculate the squared radial integral |R|^2 ---
# Based on the provided wavefunctions, the integral R = integral(r^3 * R_30 * R_31 * dr)
# can be solved analytically. The result is R = -9 * sqrt(2) * a_0 / Z.
# We need the square of its magnitude.
R_sq_val = 162 * a_0**2 / Z**2

print(f"Calculating the squared radial integral |R|^2:")
print(f"|R|^2 = 162 * a_0^2 / Z^2 = 162 * ({a_0:.4e})^2 / {Z}^2 = {R_sq_val:.4e} m^2")
print("-" * 30)

# --- Step 2: Calculate the theoretical lifetime (tau_th) ---
# The lifetime is the inverse of the Einstein A coefficient.
# The hint g2/g1=2 suggests the D2 line transition (3p_3/2 -> 3s_1/2).
# The A coefficient for this line is A_D2 = (omega^3 * e^2 * |R|^2) / (3 * pi * epsilon_0 * h_bar * c^3)
# The transition frequency is omega = 2 * pi * c / lambda.
omega = 2 * np.pi * c / lambda_wav
print(f"Calculating transition frequency omega:")
print(f"omega = 2 * pi * c / lambda = (2 * pi * {c:.4e}) / {lambda_wav:.4e} = {omega:.4e} rad/s")
print("-" * 30)

# Calculate the numerator and denominator for the lifetime formula:
# tau_th = 1 / A_D2 = (3 * pi * epsilon_0 * h_bar * c^3) / (omega^3 * e^2 * |R|^2)
numerator_tau = 3 * np.pi * epsilon_0 * h_bar * c**3
denominator_tau = omega**3 * e**2 * R_sq_val

print("Calculating the theoretical lifetime tau_th:")
print("tau_th = (3 * pi * epsilon_0 * h_bar * c^3) / (omega^3 * e^2 * |R|^2)")
print(f"tau_th = ({numerator_tau:.4e}) / (({omega:.4e})^3 * ({e:.4e})^2 * ({R_sq_val:.4e}))")
print(f"tau_th = ({numerator_tau:.4e}) / ({denominator_tau:.4e})")

tau_th = numerator_tau / denominator_tau
print(f"Theoretical lifetime tau_th = {tau_th:.4e} s = {tau_th * 1e9:.2f} ns")
print("-" * 30)

# --- Step 3: Compare theoretical and experimental lifetimes ---
ratio = tau_th / tau_exp

print("Comparing theoretical lifetime with experimental lifetime:")
print(f"Experimental lifetime tau_exp = {tau_exp * 1e9:.2f} ns")
print(f"Ratio = tau_th / tau_exp = {tau_th * 1e9:.2f} ns / {tau_exp * 1e9:.2f} ns = {ratio:.2f}")

print("\nThe theoretical lifetime is approximately 5 times the experimentally measured lifetime.")
<<<H>>>