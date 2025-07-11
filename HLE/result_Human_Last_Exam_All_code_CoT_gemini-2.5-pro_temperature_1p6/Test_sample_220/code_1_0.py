import math

def calculate_jet_speeds():
    """
    Calculates the speed of gas jets from a bursting bubble for two different sizes.
    """
    # Define physical constants
    gamma = 0.072  # Surface tension of water in N/m
    rho_air = 1.2  # Density of air in kg/m^3
    rho_water = 1000 # Density of water in kg/m^3
    g = 9.8        # Acceleration due to gravity in m/s^2

    # Define bubble diameters from the problem
    D1_mm = 2
    D2_cm = 2

    # Convert diameters to radii in meters
    R1 = (D1_mm / 1000) / 2
    R2 = (D2_cm / 100) / 2

    print("The speed 'V' of a gas jet from a bursting bubble is driven by the internal Laplace pressure.")
    print("The governing equation, from Bernoulli's principle, is V = sqrt(4 * γ / (Rc * ρ_air)), where Rc is the radius of curvature at the bubble's apex.")
    print("-" * 70)

    # --- Case 1: Small bubble (D = 2 mm) ---
    print(f"Calculating for the bubble with diameter = {D1_mm} mm (Radius R1 = {R1:.3f} m)...")
    
    # Check the regime using the Bond number
    Bo1 = (rho_water * g * R1**2) / gamma
    print(f"The Bond number Bo = (ρ_water * g * R^2) / γ = ({rho_water} * {g} * {R1:.3f}^2) / {gamma} = {Bo1:.3f}")
    print("Since Bo << 1, surface tension dominates. The bubble is hemispherical, so Rc ≈ R.")
    
    # For small bubbles, Rc = R1
    Rc1 = R1
    
    # Calculate the jet speed
    V1_sq = (4 * gamma) / (Rc1 * rho_air)
    V1 = math.sqrt(V1_sq)
    print("Jet speed V1 = sqrt(4 * γ / (R1 * ρ_air))")
    print(f"V1 = sqrt(4 * {gamma} / ({Rc1:.3f} * {rho_air})) = sqrt({V1_sq:.2f}) ≈ {V1:.2f} m/s")
    print(f"The calculated speed for the {D1_mm} mm bubble is approximately {round(V1)} m/s.")
    print("-" * 70)

    # --- Case 2: Large bubble (D = 2 cm) ---
    print(f"Calculating for the bubble with diameter = {D2_cm} cm (Radius R2 = {R2:.2f} m)...")

    # Check the regime using the Bond number
    Bo2 = (rho_water * g * R2**2) / gamma
    print(f"The Bond number Bo = (ρ_water * g * R^2) / γ = ({rho_water} * {g} * {R2:.2f}^2) / {gamma} = {Bo2:.3f}")
    print("Since Bo >> 1, gravity dominates. The bubble is flattened at the top.")
    
    # For large bubbles, Rc is approximated by the capillary length, lc
    lc_sq = gamma / (rho_water * g)
    lc = math.sqrt(lc_sq)
    print("The radius of curvature Rc is approximated by the capillary length, lc = sqrt(γ / (ρ_water * g)).")
    print(f"lc = sqrt({gamma} / ({rho_water} * {g})) = {lc:.5f} m")
    Rc2 = lc
    
    # Calculate the jet speed
    V2_sq = (4 * gamma) / (Rc2 * rho_air)
    V2 = math.sqrt(V2_sq)
    print("Jet speed V2 = sqrt(4 * γ / (lc * ρ_air))")
    print(f"V2 = sqrt(4 * {gamma} / ({Rc2:.5f} * {rho_air})) = sqrt({V2_sq:.2f}) ≈ {V2:.2f} m/s")
    print(f"The calculated speed for the {D2_cm} cm bubble is approximately {round(V2)} m/s.")
    print("-" * 70)

# Run the calculation and print the results
calculate_jet_speeds()