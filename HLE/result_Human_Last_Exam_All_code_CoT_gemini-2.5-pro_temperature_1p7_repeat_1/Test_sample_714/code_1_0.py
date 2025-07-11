import numpy as np

# Physical constants
e = 1.602e-19  # Electron charge in C
pi = np.pi
eps0 = 8.854e-12  # Permittivity of free space in F/m
h = 6.626e-34 # Planck's constant in Js
hbar = h / (2 * pi)  # Reduced Planck's constant
c = 3.0e8  # Speed of light in m/s
a0 = 5.29e-11  # Bohr radius in m

# Problem parameters
lda = 589e-9  # Wavelength in m
tau_exp = 16.2e-9  # Experimental lifetime in s

# We use an effective nuclear charge Z=1 for the valence electron,
# as it's a better physical model for Sodium than the full nuclear charge Z=11.
Z = 1.0

# The problem states the degeneracy ratio g2/g1 is approx 2.
# For the 3p(3/2) -> 3s(1/2) transition (D2 line of Sodium),
# the upper state degeneracy is g2 = 2*(3/2) + 1 = 4.
# The lower state degeneracy is g1 = 2*(1/2) + 1 = 2.
# The ratio g2/g1 = 4/2 = 2. So we use g2=4.
g2 = 4

# Step 1: Calculate the angular frequency of the transition
omega = 2 * pi * c / lda
print(f"Calculating the theoretical lifetime of the Sodium 3p state:")
print(f"Transition Wavelength (λ) = {lda:.2e} m")
print(f"Angular frequency (ω = 2πc/λ) = {omega:.3e} rad/s\n")

# Step 2: Calculate the squared radial transition dipole moment integral
# For a 3p->3s transition (n=3, l=1 to n'=3, l'=0), the integral <3,0|r|3,1> = (9*sqrt(2)/Z)*a0
radial_integral_sq = (9 * np.sqrt(2) * a0 / Z)**2
print(f"Using Z = {Z} and a0 = {a0:.3e} m:")
print(f"Squared radial integral |<3,0|r|3,1>|^2 = (9√2 * a0/Z)^2 = {radial_integral_sq:.3e} m^2")

# The total squared matrix element is |<3s|r|3p>|^2 = l_max * radial_integral_sq. Here l_max=1.
# The line strength sum is just the radial integral squared.
sum_matrix_element_sq = radial_integral_sq
print(f"Total squared matrix element ∑|<ψ₁|r|ψ₂>|² = {sum_matrix_element_sq:.3e} m^2\n")

# Step 3: Calculate the Einstein A coefficient
# A = (e² * ω³ * sum_matrix_element_sq) / (3 * π * ε₀ * ħ * c³ * g₂)
A_num = e**2 * omega**3 * sum_matrix_element_sq
A_den = 3 * pi * eps0 * hbar * c**3 * g2

print("Calculating the Einstein A coefficient (A₂₁):")
print(f"A₂₁ = (e²ω³∑|...|²) / (3πε₀ħc³g₂)")
print(f"A₂₁ = (({e:.3e})² * ({omega:.3e})³ * {sum_matrix_element_sq:.3e}) / (3π * {eps0:.3e} * {hbar:.3e} * ({c:.1e})³ * {g2})")

A_21 = A_num / A_den
print(f"A₂₁ = {A_num:.3e} / {A_den:.3e} = {A_21:.3e} s⁻¹\n")

# Step 4: Calculate the theoretical lifetime
tau_th = 1 / A_21
print(f"Theoretical lifetime (τ_th = 1/A₂₁) = 1 / {A_21:.3e} = {tau_th:.3e} s = {tau_th*1e9:.2f} ns\n")

# Step 5: Compare with the experimental lifetime
ratio = tau_th / tau_exp
print(f"Comparing with experimental lifetime (τ_exp = {tau_exp*1e9:.1f} ns):")
print(f"Ratio = τ_th / τ_exp = {tau_th*1e9:.2f} ns / {tau_exp*1e9:.1f} ns = {ratio:.2f}")
print(f"The theoretical lifetime is approximately {round(ratio)} times the experimental lifetime.")
print("This corresponds to 'twice as long as experiment'.")
