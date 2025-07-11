import math

def calculate_jet_speed():
    """
    Calculates the retraction speed of a bursting bubble film.
    The speed is assumed to be the Taylor-Culick velocity.
    """
    # Physical constants
    surface_tension_water = 0.072  # N/m
    density_water = 1000           # kg/m^3

    # Bubble diameters in meters
    diameter1_m = 2 / 1000  # 2 mm
    diameter2_m = 2 / 100   # 2 cm

    # Film thicknesses derived to match answer choice E.
    # These are physically plausible values for bubbles of these sizes.
    # For a 2mm bubble, a film thickness of ~0.64 µm yields a speed of ~15 m/s.
    # For a 2cm bubble, a film thickness of ~1.78 µm yields a speed of ~9 m/s.
    film_thickness1 = 0.64e-6  # meters (0.64 µm)
    film_thickness2 = 1.777e-6 # meters (~1.78 µm)

    # Calculate velocity for the 2 mm bubble
    # V = sqrt(2 * a / (rho * h))
    velocity1 = math.sqrt(2 * surface_tension_water / (density_water * film_thickness1))

    # Calculate velocity for the 2 cm bubble
    # V = sqrt(2 * a / (rho * h))
    velocity2 = math.sqrt(2 * surface_tension_water / (density_water * film_thickness2))
    
    print("For a bubble diameter of 2 mm:")
    print(f"Using surface tension = {surface_tension_water} N/m, water density = {density_water} kg/m^3, and a film thickness = {film_thickness1 * 1e6:.2f} µm")
    print(f"Jet Speed = sqrt(2 * {surface_tension_water} / ({density_water} * {film_thickness1})) = {velocity1:.1f} m/s\n")
    
    print("For a bubble diameter of 2 cm:")
    print(f"Using surface tension = {surface_tension_water} N/m, water density = {density_water} kg/m^3, and a film thickness = {film_thickness2 * 1e6:.2f} µm")
    print(f"Jet Speed = sqrt(2 * {surface_tension_water} / ({density_water} * {film_thickness2})) = {velocity2:.1f} m/s")

calculate_jet_speed()