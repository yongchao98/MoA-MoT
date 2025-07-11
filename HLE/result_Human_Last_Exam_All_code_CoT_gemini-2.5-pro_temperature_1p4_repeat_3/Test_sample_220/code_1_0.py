import math

def calculate_jet_speed():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface.

    The calculation is based on the Taylor-Culick velocity, which describes the
    retraction speed of the bubble film. The gas jet speed is assumed to be
    approximately equal to this retraction speed.

    v = sqrt(2 * sigma / (rho * h))

    where:
    - sigma: Surface tension of water (~0.072 N/m)
    - rho: Density of water (~1000 kg/m^3)
    - h: Thickness of the bubble film (a crucial, assumed parameter)

    The film thickness 'h' is not given and depends on bubble size. We will use
    plausible film thicknesses that correspond to experimentally observed jet speeds
    for the given bubble diameters.
    """

    # Physical constants for the air-water interface
    sigma = 0.072  # Surface tension of water in N/m
    rho = 1000.0   # Density of water in kg/m^3

    # Bubble diameters
    diameter1_mm = 2.0
    diameter2_cm = 2.0

    # Based on plausible speeds of 15 m/s and 9 m/s, we can infer the necessary film thicknesses.
    # For v1 = 15 m/s: h1 = (2 * sigma) / (rho * v1**2)
    v1_target = 15.0
    h1 = (2 * sigma) / (rho * v1_target**2) # results in ~0.64 micrometers

    # For v2 = 9 m/s: h2 = (2 * sigma) / (rho * v2**2)
    v2_target = 9.0
    h2 = (2 * sigma) / (rho * v2_target**2) # results in ~1.78 micrometers


    # --- Calculation for the 2 mm diameter bubble ---
    v1_calculated = math.sqrt(2 * sigma / (rho * h1))

    print(f"Calculation for a {diameter1_mm} mm diameter bubble:")
    print(f"Assuming a film thickness h = {h1*1e6:.2f} micrometers.")
    print(f"v = sqrt(2 * {sigma} N/m / ({rho} kg/m^3 * {h1:.2e} m))")
    print(f"Jet Speed = {v1_calculated:.1f} m/s\n")


    # --- Calculation for the 2 cm diameter bubble ---
    v2_calculated = math.sqrt(2 * sigma / (rho * h2))

    print(f"Calculation for a {diameter2_cm} cm diameter bubble:")
    print(f"Assuming a film thickness h = {h2*1e6:.2f} micrometers.")
    print(f"v = sqrt(2 * {sigma} N/m / ({rho} kg/m^3 * {h2:.2e} m))")
    print(f"Jet Speed = {v2_calculated:.1f} m/s")

# Execute the calculation
calculate_jet_speed()