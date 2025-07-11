import math

# --- Step 1: Define Constants and Parameters ---
q = 1.602e-19  # Elementary charge in C
eps0 = 8.854e-12 # Permittivity of free space in F/m
h_bar = 1.05457e-34 # Reduced Planck constant in J*s
c = 2.998e8   # Speed of light in m/s
a0 = 5.29e-11   # Bohr radius in m
pi = math.pi

# Problem specific parameters
Z = 11.0  # Nuclear charge
lambda_ = 589e-9  # Wavelength in m
tau_exp_ns = 16.2 # Experimental lifetime in ns
tau_exp = tau_exp_ns * 1e-9 # Experimental lifetime in s
# Per the problem statement, we use the effective degeneracy for the upper state.
# Orbital degeneracy would be 2*l+1 = 3. Hint suggests using g2/g1=2. Since g1=1, g2=2.
g2 = 2
l_greater = 1 # For a p->s transition (l=1 to l=0), the greater l is 1

# --- Step 2: Calculate the Radial Integral Squared ---
# The radial integral R_integral = integral(r^3 * R_3,0 * R_3,1 dr)
# Standard calculation for hydrogenic wavefunctions yields: R_integral = -9 * sqrt(2) * a0 / Z
R_integral = -9 * math.sqrt(2) * a0 / Z
R_integral_sq = R_integral**2

# --- Step 3: Calculate the Sum of Squared Matrix Elements ---
# Sum |<psi_1|r|psi_2>|^2 = l_greater * R_integral_sq
sum_matrix_element_sq = l_greater * R_integral_sq

# --- Step 4: Calculate the Einstein A Coefficient (Spontaneous Emission Rate) ---
# A21 = (8 * pi^2 * q^2) / (3 * eps0 * h_bar * lambda^3) * (1/g2) * sum_matrix_element_sq
term1 = (8 * pi**2 * q**2) / (3 * eps0 * h_bar * lambda_**3)
A21 = term1 * (1 / g2) * sum_matrix_element_sq

# Alternatively, using the derived final expression for tau_th
num = 3 * eps0 * h_bar * lambda_**3 * Z**2
den = 648 * pi**2 * q**2 * a0**2
tau_th_calc = num / den

# --- Step 5: Calculate Theoretical Lifetime and Compare ---
tau_th = 1 / A21
ratio = tau_th / tau_exp

# --- Step 6: Print the Calculation and Result ---
print("Calculating the theoretical lifetime (tau_th) of the Sodium-23 3p state.")
print("\nThe formula for the lifetime is: tau_th = (3 * eps0 * h_bar * lambda^3 * Z^2) / (648 * pi^2 * q^2 * a0^2)")
print("\nPlugging in the values:")
print(f"tau_th = (3 * {eps0:.3e} * {h_bar:.3e} * ({lambda_:.3e})^3 * {Z:.0f}^2) / (648 * {pi:.4f}^2 * ({q:.3e})^2 * ({a0:.3e})^2)")

numerator = 3 * eps0 * h_bar * (lambda_**3) * (Z**2)
denominator = 648 * (pi**2) * (q**2) * (a0**2)

print(f"\nNumerator = {numerator:.4e}")
print(f"Denominator = {denominator:.4e}")

tau_th_ns = tau_th * 1e9
print(f"\nTheoretical lifetime (tau_th) = {numerator:.4e} / {denominator:.4e} = {tau_th:.4e} s")
print(f"Theoretical lifetime (tau_th) = {tau_th_ns:.1f} ns")

print(f"\nThe experimentally measured lifetime (tau_exp) is {tau_exp_ns:.1f} ns.")
print(f"The ratio of theoretical to experimental lifetime is {tau_th_ns:.1f} ns / {tau_exp_ns:.1f} ns = {ratio:.2f}")

print("\nThis ratio is approximately 10.")
print("The theoretical lifetime is about ten times as long as the experimental result.")
<<<G>>>