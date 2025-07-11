import math

def calculate_jet_speed():
    """
    Calculates the speed of gas jets from bursting bubbles of two different sizes.

    The model uses Bernoulli's principle driven by Laplace pressure. It accounts
    for the different pressure regimes for small (capillary-dominated) and
    large (gravity-dominated) bubbles.
    """
    # Physical constants
    sigma = 0.072  # Surface tension of water (N/m)
    rho_water = 1000  # Density of water (kg/m^3)
    rho_air = 1.2  # Density of air (kg/m^3)
    g = 9.8  # Acceleration due to gravity (m/s^2)

    # Diameters of the bubbles
    d1 = 0.002  # 2 mm in meters
    d2 = 0.02   # 2 cm in meters
    
    # Radii of the bubbles
    R1 = d1 / 2
    R2 = d2 / 2

    print("--- Model Parameters ---")
    print(f"Surface tension of water (sigma): {sigma} N/m")
    print(f"Density of air (rho_air): {rho_air} kg/m^3")
    print(f"Density of water (rho_water): {rho_water} kg/m^3")

    # 1. Calculate the capillary length
    l_c = math.sqrt(sigma / (rho_water * g))
    print(f"\nCapillary length (l_c) = sqrt({sigma} / ({rho_water} * {g})) = {l_c:.4f} m or {l_c*1000:.2f} mm")
    print("This length separates small (R < l_c) from large (R > l_c) bubbles.")

    # 2. Calculation for the 2 mm diameter bubble
    print("\n--- Case 1: Bubble Diameter = 2 mm ---")
    print(f"Radius R1 = {R1*1000:.1f} mm. Since R1 < l_c, bubble is capillary-dominated.")
    
    # Laplace pressure for a small hemispherical bubble
    delta_P1 = (2 * sigma) / R1
    print("Equation for pressure: delta_P1 = 2 * sigma / R1")
    print(f"delta_P1 = 2 * {sigma} / {R1} = {delta_P1:.2f} Pa")
    
    # Jet speed from Bernoulli's principle
    v1 = math.sqrt((2 * delta_P1) / rho_air)
    print("Equation for speed: v1 = sqrt(2 * delta_P1 / rho_air)")
    print(f"v1 = sqrt(2 * {delta_P1:.2f} / {rho_air}) = {v1:.1f} m/s")

    # 3. Calculation for the 2 cm diameter bubble
    print("\n--- Case 2: Bubble Diameter = 2 cm ---")
    print(f"Radius R2 = {R2*1000:.1f} mm. Since R2 > l_c, bubble is gravity-dominated.")
    
    # Laplace pressure for a large flattened bubble, curvature ~ 1/l_c
    delta_P2 = (2 * sigma) / l_c
    print("Equation for pressure: delta_P2 = 2 * sigma / l_c")
    print(f"delta_P2 = 2 * {sigma} / {l_c:.4f} = {delta_P2:.2f} Pa")
    
    # Jet speed from Bernoulli's principle
    v2 = math.sqrt((2 * delta_P2) / rho_air)
    print("Equation for speed: v2 = sqrt(2 * delta_P2 / rho_air)")
    print(f"v2 = sqrt(2 * {delta_P2:.2f} / {rho_air}) = {v2:.1f} m/s")
    
    print("\n--- Final Results ---")
    print(f"The calculated speed for a 2 mm bubble is approximately {round(v1)} m/s.")
    print(f"The calculated speed for a 2 cm bubble is approximately {round(v2)} m/s.")

calculate_jet_speed()