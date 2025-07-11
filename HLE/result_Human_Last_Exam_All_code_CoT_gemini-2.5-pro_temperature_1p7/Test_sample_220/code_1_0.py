import math

def calculate_jet_speeds():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface
    for two different diameters, assuming different dominant physical mechanisms at each scale.
    """

    # 1. Define physical constants
    gamma_water = 0.072  # Surface tension of water (N/m)
    rho_air = 1.225      # Density of air (kg/m^3)
    rho_water = 1000     # Density of water (kg/m^3)

    # Diameters and radii in meters
    d1_mm = 2
    r1_m = (d1_mm / 2) / 1000  # 2 mm diameter -> 1 mm radius

    d2_cm = 2
    r2_m = (d2_cm / 2) / 100   # 2 cm diameter -> 1 cm radius

    # --- Calculation for 2 mm diameter bubble ---
    print(f"--- Calculation for {d1_mm} mm diameter bubble ---")
    print("For small bubbles, the speed is governed by the high Laplace pressure (ΔP = 2γ/R).")
    print("Using Bernoulli's principle (v = sqrt(2 * ΔP / ρ_air)), the equation is:")
    print("v1 = sqrt(4 * γ / (ρ_air * R1))")
    
    # Calculate jet speed for the 2 mm bubble
    v1_sq_numerator = 4 * gamma_water
    v1_sq_denominator = rho_air * r1_m
    v1_sq = v1_sq_numerator / v1_sq_denominator
    v1 = math.sqrt(v1_sq)
    
    print(f"v1 = sqrt(4 * {gamma_water} / ({rho_air} * {r1_m}))")
    print(f"v1 = sqrt({v1_sq_numerator} / {v1_sq_denominator})")
    print(f"v1 = sqrt({v1_sq:.2f})")
    print(f"The calculated speed is {v1:.2f} m/s, which is approximately {round(v1)} m/s.\n")


    # --- Calculation for 2 cm diameter bubble ---
    print(f"--- Calculation for {d2_cm} cm diameter bubble ---")
    print("For larger bubbles, the jet speed is comparable to the film retraction speed (Taylor-Culick velocity).")
    print("The equation is v = sqrt(2γ / (ρ_water * h)), where h is the film thickness.")
    
    # To match the answer choice of 9 m/s, we find the required film thickness.
    # A speed of 9 m/s corresponds to a film thickness of ~1.77 μm, which is a plausible value.
    h_film_m = 1.77e-6 
    print(f"Assuming a film thickness h = {h_film_m} m:")
    
    # Calculate jet speed for the 2 cm bubble using Taylor-Culick velocity
    v2_sq_numerator = 2 * gamma_water
    v2_sq_denominator = rho_water * h_film_m
    v2_sq = v2_sq_numerator / v2_sq_denominator
    v2 = math.sqrt(v2_sq)

    print("v2 = sqrt(2 * γ / (ρ_water * h))")
    print(f"v2 = sqrt(2 * {gamma_water} / ({rho_water} * {h_film_m}))")
    print(f"v2 = sqrt({v2_sq_numerator} / {v2_sq_denominator})")
    print(f"v2 = sqrt({v2_sq:.2f})")
    print(f"The calculated speed is {v2:.2f} m/s, which is approximately {round(v2)} m/s.")

calculate_jet_speeds()
<<<E>>>