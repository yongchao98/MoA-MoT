import numpy as np

# Step 1: Define physical constants and given parameters
e = 1.602e-19  # Elementary charge in C
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m
h_bar = 1.054e-34   # Reduced Planck constant in J*s
c = 3.0e8         # Speed of light in m/s
a_0 = 5.29e-11    # Bohr radius in m

# Given experimental data
lambda_nm = 589      # Wavelength in nm
lambda_m = lambda_nm * 1e-9 # Wavelength in m
tau_exp_ns = 16.2    # Experimental lifetime in ns
tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in s

# Step 2: Set model parameters
# For a Sodium valence electron, the core electrons screen the nucleus.
# We model it as a hydrogenic atom with an effective nuclear charge Z_eff = 1.
Z_eff = 1
# The degeneracy of the upper state (3p, l=1) is g2 = 2*l + 1
g2 = 3

print("--- Calculation Plan ---")
print("1. Calculate angular frequency (omega) from the wavelength.")
print("2. Calculate the squared radial integral |<3s|r|3p>|^2 using the hydrogenic model with Z_eff = 1.")
print("3. Calculate the Einstein A coefficient (A21).")
print("4. Calculate the theoretical lifetime (tau_th).")
print("5. Compare theoretical and experimental lifetimes.")
print("-" * 25)

# Step 3: Calculate angular frequency omega
omega = 2 * np.pi * c / lambda_m
print(f"Step 1: Angular Frequency")
print(f"ω = 2 * π * c / λ = 2 * π * {c:.2e} / {lambda_m:.2e} = {omega:.3e} rad/s")

# Step 4: Calculate the squared radial integral
# For a hydrogenic n,l -> n,l-1 transition, the integral <R_nl-1|r|R_nl> = (3n/2Z) * sqrt(n^2 - l^2) * a0
# For 3p -> 3s, n=3, l=1.
radial_integral = (3 * 3 / (2 * Z_eff)) * np.sqrt(3**2 - 1**2) * a_0
radial_integral_sq = radial_integral**2
# The result is equivalent to 162 * (a_0/Z_eff)^2
check_integral_sq = 162 * (a_0 / Z_eff)**2

print(f"\nStep 2: Squared Radial Integral")
print(f"I = <R_3s|r|R_3p> = (3*3/(2*{Z_eff})) * sqrt(3^2 - 1^2) * a_0 = {radial_integral:.3e} m")
print(f"|I|^2 = ({radial_integral:.3e})^2 = {radial_integral_sq:.3e} m^2")


# Step 5: Calculate the Einstein A coefficient (A21)
# A21 = (e^2 * omega^3) / (3 * pi * epsilon_0 * h_bar * c^3 * g2) * |I|^2
# Note: The formula in many textbooks simplifies to (4 * alpha * omega^3) / (3 * c^2 * g2) * |I|^2 where alpha is the fine-structure constant.
# We will use the fundamental constants directly.
numerator = e**2 * omega**3 * radial_integral_sq
denominator = 3 * np.pi * epsilon_0 * h_bar * c**3 * g2
A21 = numerator / denominator

print(f"\nStep 3: Einstein A Coefficient (A21)")
print(f"A21 = (e^2 * ω^3 * |I|^2) / (3 * π * ε₀ * ħ * c^3 * g₂) ")
print(f"A21 = (({e:.3e})^2 * ({omega:.3e})^3 * {radial_integral_sq:.3e}) / (3 * π * {epsilon_0:.3e} * {h_bar:.3e} * ({c:.2e})^3 * {g2})")
print(f"A21 = {A21:.3e} s^-1")


# Step 6: Calculate theoretical lifetime
tau_th_s = 1 / A21
tau_th_ns = tau_th_s * 1e9

print(f"\nStep 4: Theoretical Lifetime (τ_th)")
print(f"τ_th = 1 / A21 = 1 / {A21:.3e} = {tau_th_s:.3e} s")
print(f"τ_th = {tau_th_ns:.3f} ns")

# Step 7: Compare and conclude
ratio = tau_th_s / tau_exp_s
print(f"\nStep 5: Comparison")
print(f"Experimental Lifetime τ_exp = {tau_exp_ns} ns")
print(f"Theoretical Lifetime τ_th = {tau_th_ns:.3f} ns")
print(f"The ratio τ_th / τ_exp = {tau_th_ns:.3f} / {tau_exp_ns} = {ratio:.3f}")

print("\nThe theoretical lifetime is approximately 0.1 times the experimental lifetime.")
print("This corresponds to answer choice C.")
<<<C>>>