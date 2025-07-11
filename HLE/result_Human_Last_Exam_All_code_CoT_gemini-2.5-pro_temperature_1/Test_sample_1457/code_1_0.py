import math

# This script calculates the X-ray transmission through the detector's Be window
# to find the lightest element that can be detected.

# --- Physical Constants and Parameters ---
# Density of Beryllium (ρ) in g/cm^3
rho_be = 1.85
# Thickness of the Be window (x) in cm (100 µm = 0.01 cm)
x_be = 0.01

# --- Element Data ---
# We will check the K-alpha X-ray line for each element.
# Data includes the element's name, K-alpha energy (keV), and the
# mass attenuation coefficient (μ) of Be at that energy (cm^2/g).
elements_data = [
    {
        "name": "Sodium (Na)",
        "energy_keV": 1.041,
        "mu_be_cm2_g": 522
    },
    {
        "name": "Magnesium (Mg)",
        "energy_keV": 1.254,
        "mu_be_cm2_g": 341
    },
    {
        "name": "Silicon (Si)",
        "energy_keV": 1.740,
        "mu_be_cm2_g": 167
    },
    {
        "name": "Calcium (Ca)",
        "energy_keV": 3.692,
        "mu_be_cm2_g": 27.2
    }
]

print("Calculating X-ray transmission through a 100 µm Beryllium window.")
print("The lightest element with a transmission greater than ~1% is considered detectable.\n")

# Loop through each element to calculate its transmission
for element in elements_data:
    name = element["name"]
    mu = element["mu_be_cm2_g"]

    # Calculate transmission using the Beer-Lambert law
    transmission_fraction = math.exp(-mu * rho_be * x_be)
    transmission_percent = transmission_fraction * 100

    print(f"--- Analysis for {name} ---")
    print(f"The equation for transmission is: T(%) = 100 * exp(-μ * ρ * x)")
    # Printing the equation with the specific numbers for this element
    print(f"Calculation: T(%) = 100 * exp(-{mu} * {rho_be} * {x_be})")
    print(f"Result: The transmission is {transmission_percent:.4f}%\n")

print("Conclusion: Based on the calculations, Silicon (Si) is the lightest element")
print("from the list with a significant transmission percentage, making it reliably detectable.")
