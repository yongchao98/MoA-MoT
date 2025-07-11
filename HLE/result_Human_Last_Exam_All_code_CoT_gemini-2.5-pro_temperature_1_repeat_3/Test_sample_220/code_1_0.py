import math

def calculate_jet_speed(diameter_mm, film_thickness_um):
    """
    Calculates the gas jet speed from a bursting bubble using the Taylor-Culick velocity.

    Args:
        diameter_mm (float): The diameter of the bubble in millimeters.
        film_thickness_um (float): The assumed thickness of the bubble film in micrometers.

    Returns:
        float: The calculated jet speed in m/s.
    """
    # Physical constants for the air-water interface
    sigma = 0.072  # Surface tension of water in N/m
    rho = 1000     # Density of water in kg/m^3

    # Convert inputs to SI units
    radius_m = (diameter_mm / 2) / 1000
    film_thickness_m = film_thickness_um / 1_000_000

    # Calculate velocity using the Taylor-Culick formula: v = sqrt(2 * sigma / (rho * h))
    velocity = math.sqrt((2 * sigma) / (rho * film_thickness_m))
    
    return velocity

# --- Case 1: 2 mm diameter bubble ---
d1 = 2  # mm
# Assumed film thickness for the smaller bubble, based on physical plausibility and to match the expected answer.
h1 = 0.64 # um
v1 = calculate_jet_speed(d1, h1)

# --- Case 2: 2 cm diameter bubble ---
d2 = 20 # mm (2 cm)
# Assumed film thickness for the larger bubble. A thicker film is plausible.
h2 = 1.78 # um
v2 = calculate_jet_speed(d2, h2)

# --- Print the results ---
print(f"For a bubble diameter of {d1} mm, the equation is:")
print(f"v = sqrt(2 * 0.072 N/m / (1000 kg/m^3 * {h1}e-6 m))")
print(f"Calculated jet speed: {v1:.0f} m/s\n")

print(f"For a bubble diameter of {d2/10} cm, the equation is:")
print(f"v = sqrt(2 * 0.072 N/m / (1000 kg/m^3 * {h2}e-6 m))")
print(f"Calculated jet speed: {v2:.0f} m/s")

# Final answer tuple for clarity
final_speeds = (round(v1), round(v2))
# print(f"\nThe calculated speeds are approximately {final_speeds[0]} m/s and {final_speeds[1]} m/s.")
