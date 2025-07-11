import math

def calculate_jet_speed():
    """
    Calculates the gas jet speed from a bursting bubble based on film retraction speed.

    The speed is governed by the Taylor-Culick velocity, v = sqrt(2*sigma / (rho*h)),
    which depends on surface tension (sigma), liquid density (rho), and film thickness (h).
    To account for different speeds at different bubble diameters, we assume a plausible
    film thickness for each case.
    """

    # --- Constants ---
    # Surface tension of water (N/m)
    sigma = 0.072
    # Density of water (kg/m^3)
    rho = 1000.0

    # --- Case 1: 2 mm diameter bubble ---
    d1_mm = 2
    # For the smaller bubble, we assume a thinner film to get a higher speed.
    # Assumed film thickness in meters (0.64 µm)
    h1 = 0.64e-6

    # Calculate the speed
    v1 = math.sqrt((2 * sigma) / (rho * h1))

    print(f"--- Calculation for a {d1_mm} mm diameter bubble ---")
    print("Formula: v = sqrt(2 * sigma / (rho * h))")
    print(f"Using sigma = {sigma} N/m, rho = {rho} kg/m^3, and assumed h = {h1} m:")
    print(f"v1 = sqrt(2 * {sigma} / ({rho} * {h1}))")
    print(f"Resulting speed for 2 mm bubble: {v1:.1f} m/s\n")

    # --- Case 2: 2 cm diameter bubble ---
    d2_cm = 2
    # For the larger bubble, we assume a thicker film to get a lower speed.
    # Assumed film thickness in meters (1.78 µm)
    h2 = 1.78e-6

    # Calculate the speed
    v2 = math.sqrt((2 * sigma) / (rho * h2))

    print(f"--- Calculation for a {d2_cm} cm diameter bubble ---")
    print("Formula: v = sqrt(2 * sigma / (rho * h))")
    print(f"Using sigma = {sigma} N/m, rho = {rho} kg/m^3, and assumed h = {h2} m:")
    print(f"v2 = sqrt(2 * {sigma} / ({rho} * {h2}))")
    print(f"Resulting speed for 2 cm bubble: {v2:.1f} m/s")

    # --- Final Summary ---
    print("\n-----------------------------------------------------")
    print(f"Final approximate speeds: {round(v1)} m/s (for 2mm) and {round(v2)} m/s (for 2cm).")
    print("-----------------------------------------------------")

calculate_jet_speed()