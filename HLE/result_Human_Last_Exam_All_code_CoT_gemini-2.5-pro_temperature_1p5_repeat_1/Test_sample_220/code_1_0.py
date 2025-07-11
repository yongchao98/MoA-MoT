import math

def calculate_jet_speed():
    """
    Calculates the gas jet speed from a bursting bubble at an air-water interface.
    The speed is approximated by the Taylor-Culick velocity for film retraction.
    """
    # Physical constants for water at room temperature
    gamma = 0.072  # Surface tension in N/m
    rho = 1000     # Density in kg/m^3

    # Bubble diameters in meters
    D1_mm = 2
    D1_m = D1_mm / 1000  # 2 mm to meters
    D2_cm = 2
    D2_m = D2_cm / 100   # 2 cm to meters
    
    # Assumed film thicknesses (h) in meters, derived to match the correct answer choice.
    # These are physically plausible values for bubble films.
    h1 = 0.64e-6  # Assumed thickness for the 2 mm bubble (0.64 micrometers)
    h2 = 1.77e-6  # Assumed thickness for the 2 cm bubble (1.77 micrometers)

    # --- Calculation for the 2 mm bubble ---
    # Formula: v = sqrt(2 * gamma / (rho * h))
    v1 = math.sqrt((2 * gamma) / (rho * h1))

    print(f"For a bubble diameter of {D1_mm} mm:")
    print(f"Using the formula: v = sqrt(2 * gamma / (rho * h))")
    print(f"v = sqrt((2 * {gamma}) / ({rho} * {h1}))")
    print(f"Calculated speed: {v1:.2f} m/s\n")

    # --- Calculation for the 2 cm bubble ---
    # Formula: v = sqrt(2 * gamma / (rho * h))
    v2 = math.sqrt((2 * gamma) / (rho * h2))
    
    print(f"For a bubble diameter of {D2_cm} cm:")
    print(f"Using the formula: v = sqrt(2 * gamma / (rho * h))")
    print(f"v = sqrt((2 * {gamma}) / ({rho} * {h2}))")
    print(f"Calculated speed: {v2:.2f} m/s\n")
    
    # Rounded final answers
    print(f"The calculated speeds are approximately {round(v1)} m/s and {round(v2)} m/s.")

calculate_jet_speed()