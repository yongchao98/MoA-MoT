import cmath

# In the Feynman path integral formalism, the bare Green's function G_0
# for a non-interacting particle with energy eigenvalue epsilon_k
# is expressed in the frequency-momentum domain.

# --- Parameters ---
# Let's choose some example values for the variables.
# The units would typically be electron-volts (eV) in solid-state physics.

# The energy of the single-particle state (e.g., for a specific momentum k)
epsilon_k = 1.5

# The frequency (probing energy) at which we evaluate the Green's function
omega = 1.4

# A small positive infinitesimal value for the imaginary part
eta = 0.1

# --- Functional Dependence ---
# The dependence of the bare Green's function G_0 on the energy
# eigenvalue epsilon_k is given by the following formula for the
# retarded Green's function:
#
# G_0(k, omega) = 1 / (omega - epsilon_k + i*eta)

print("The functional form of the bare Green's function G_0 is:")
print("  G_0(k, ω) = 1 / (ω - ε_k + iη)")
print("-" * 40)
print("Substituting the example numerical values:")
print(f"  Probing Frequency (ω) = {omega}")
print(f"  Single-particle Energy (ε_k) = {epsilon_k}")
print(f"  Infinitesimal (η) = {eta}")
print("-" * 40)

# We construct the final equation string using the values.
# This clearly shows each number in its place within the equation.
print("The final equation is:")
print(f"  G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")

# Let's also compute the numerical result for completeness.
denominator = (omega - epsilon_k) + 1j * eta
g0_value = 1 / denominator

print("\nCalculating the result:")
print(f"  G_0 = 1 / ({denominator.real} + {denominator.imag}i)")
print(f"  G_0 ≈ {g0_value.real:.2f} {g0_value.imag:+.2f}i")
