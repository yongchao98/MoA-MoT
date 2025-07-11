import math

# --- Define the numerical coefficients from the derivation ---
# The sum of squares of orientational factors (1^2 + 1^2 + (-2)^2) gives 6.
# The energy denominator (-2*h_bar*omega_0) gives a factor of -1/2.
# The square of the interaction prefactor (1/(4*pi))^2 gives 1/(16*pi^2).
# The square of the matrix element prefactor (1/2)^2 gives 1/4.
# Total numeric coefficient: 6 * (-1/2) * (1/16) * (1/4) would be incorrect this way.
# Let's use the final derived numbers.

# Final coefficient in the numerator
numerator_coefficient = 3

# Final coefficient in the denominator
denominator_coefficient = 64

# --- Print the explanation and the final formula ---
print("The ground state energy shift (ΔE₀) due to the Coulomb interaction is found using second-order perturbation theory.")
print("This interaction is also known as the Van der Waals or London dispersion force.")
print("The final expression for the energy shift is:")
print() # Add a blank line for spacing

# Construct and print the formatted equation
# The equation is: - (3 * e^4 * h_bar) / (64 * pi^2 * m^2 * omega_0^3 * R^6)

# Numerator string
numerator_str = f"  {numerator_coefficient} \u00B7 e\u2074 \u00B7 \u0127"

# Denominator string
denominator_str = f"{denominator_coefficient} \u00B7 \u03c0\u00b2 \u00B7 m\u00b2 \u00B7 \u03c9\u2080\u00b3 \u00B7 R\u2076"

# Create the fraction bar
bar_length = max(len(numerator_str), len(denominator_str))
fraction_bar = "\u2014" * bar_length

# Print the equation part by part
print(f"ΔE₀ = - {numerator_str}")
print(f"        {fraction_bar}")
print(f"        {denominator_str}")
print()
print("Where:")
print("e       = electron charge")
print("\u0127       = reduced Planck constant")
print("\u03c0       = pi")
print("m       = mass of the oscillator")
print("\u03c9\u2080      = angular frequency of the oscillator")
print("R       = distance between the oscillators")