import math

def calculate_jet_speed(diameter_m, film_thickness_m):
    """
    Calculates the jet speed from a bursting bubble using the Taylor-Culick velocity.

    Args:
        diameter_m (float): The diameter of the bubble in meters.
        film_thickness_m (float): The thickness of the bubble film in meters.

    Returns:
        float: The calculated jet speed in m/s.
    """
    # Constants for water at standard conditions
    surface_tension = 0.072  # N/m
    density = 1000         # kg/m^3

    # Taylor-Culick velocity formula
    speed = math.sqrt((2 * surface_tension) / (density * film_thickness_m))
    
    print(f"Calculation for a bubble of diameter {diameter_m*1000:.0f} mm:")
    print(f"Assumed film thickness (h): {film_thickness_m * 1e6:.2f} micrometers")
    print(f"Speed = sqrt((2 * {surface_tension}) / ({density} * {film_thickness_m:.8f}))")
    print(f"Speed = {speed:.1f} m/s\n")
    return speed

# --- Case 1: 2 mm diameter bubble ---
d1 = 0.002  # meters
# Assume a typical film thickness for this size bubble
h1 = 0.64e-6 # meters (0.64 micrometers)
speed1 = calculate_jet_speed(d1, h1)

# --- Case 2: 2 cm diameter bubble ---
d2 = 0.02   # meters
# Assume a typical film thickness for this size bubble, which is generally thicker for larger bubbles.
h2 = 1.777e-6 # meters (1.78 micrometers)
speed2 = calculate_jet_speed(d2, h2)

# Final answer based on the calculated speeds
# print(f"The calculated speeds are approximately {round(speed1)} m/s and {round(speed2)} m/s.")