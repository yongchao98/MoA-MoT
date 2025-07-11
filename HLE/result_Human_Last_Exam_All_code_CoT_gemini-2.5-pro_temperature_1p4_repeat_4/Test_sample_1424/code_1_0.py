import cmath

# 1. Define the parameters for the calculation.
# We'll use example values for the frequency (omega), the single-particle energy (epsilon_k),
# and the infinitesimal broadening term (eta).
omega = 2.5  # Energy of the probe, in some arbitrary units
epsilon_k = 2.0  # Single-particle energy eigenvalue
eta = 0.1  # A small positive constant for causality

# 2. State the formula being used.
print("The retarded bare Green's function G_0 is a function of frequency ω and the single-particle energy ϵ_k:")
print("G_0(ω, ϵ_k) = 1 / (ω - ϵ_k + iη)\n")

# 3. Substitute the numerical values into the formula and print the step-by-step calculation.
print("For our chosen values:")
print(f"ω = {omega}")
print(f"ϵ_k = {epsilon_k}")
print(f"η = {eta}\n")

print("The equation becomes:")
# Python uses 'j' for the imaginary unit, but we print 'i' for convention.
print(f"G_0 = 1 / (({omega} - {epsilon_k}) + {eta}i)")

# 4. Perform the calculation.
# The denominator is a complex number.
# Python uses 'j' for the imaginary unit.
denominator_real = omega - epsilon_k
denominator_imag = eta
denominator = complex(denominator_real, denominator_imag)

# Calculate G_0
G0 = 1 / denominator

# 5. Print the intermediate step and the final result.
print(f"G_0 = 1 / ({denominator.real} + {denominator.imag}i)")
print(f"\nFinal Result:")
print(f"G_0 = {G0}")
print(f"G_0 ≈ {G0.real:.4f} + ({G0.imag:.4f})i")