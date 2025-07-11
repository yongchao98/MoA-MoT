import math

def calculate_jet_speed():
    """
    Calculates the speed of gas jets from bursting bubbles of two different sizes
    at an air-water interface, considering the dominant physical regime for each size.
    """
    # Physical constants
    sigma = 0.072  # Surface tension of water (N/m)
    rho_water = 1000  # Density of water (kg/m^3)
    rho_air = 1.2  # Density of air (kg/m^3)
    g = 9.8  # Acceleration due to gravity (m/s^2)

    # Bubble diameters and radii
    D1 = 0.002  # 2 mm in meters
    R1 = D1 / 2
    D2 = 0.02  # 2 cm in meters
    R2 = D2 / 2

    # --- Step 1: Determine the physical regimes by calculating the capillary length ---
    capillary_length = math.sqrt(sigma / (rho_water * g))
    print(f"Physical Constants:")
    print(f"  Surface Tension of Water (σ): {sigma} N/m")
    print(f"  Density of Air (ρ_air): {rho_air} kg/m^3")
    print(f"  Density of Water (ρ_water): {rho_water} kg/m^3")
    print(f"  Gravity (g): {g} m/s^2\n")

    print(f"First, we calculate the capillary length to determine which forces (surface tension or gravity) are dominant.")
    print(f"Capillary Length (λ_c) = sqrt(σ / (ρ_water * g)) = sqrt({sigma} / ({rho_water} * {g})) = {capillary_length:.4f} m or {capillary_length*1000:.1f} mm\n")

    # --- Step 2: Calculate jet speed for the 2 mm bubble ---
    print("--- Case 1: Bubble Diameter = 2 mm ---")
    print(f"Bubble radius R1 = {R1*1000:.1f} mm. Since R1 < λ_c ({R1*1000:.1f} mm < {capillary_length*1000:.1f} mm), surface tension dominates.")
    
    # Laplace pressure for a small bubble
    delta_p1 = (2 * sigma) / R1
    print(f"The excess pressure ΔP1 is given by the Laplace equation: ΔP1 = 2 * σ / R1")
    print(f"ΔP1 = 2 * {sigma} / {R1} = {delta_p1:.2f} Pa")
    
    # Jet speed from Bernoulli's principle
    v_jet1 = math.sqrt((2 * delta_p1) / rho_air)
    print(f"The jet speed V1 is calculated using Bernoulli's principle: V1 = sqrt(2 * ΔP1 / ρ_air)")
    print(f"V1 = sqrt(2 * {delta_p1:.2f} / {rho_air}) = {v_jet1:.2f} m/s\n")
    
    # --- Step 3: Calculate jet speed for the 2 cm bubble ---
    print("--- Case 2: Bubble Diameter = 2 cm ---")
    print(f"Bubble radius R2 = {R2*1000:.1f} mm. Since R2 > λ_c ({R2*1000:.1f} mm > {capillary_length*1000:.1f} mm), gravity flattens the bubble.")

    # Pressure for a large, gravity-flattened bubble
    delta_p2 = 2 * math.sqrt(sigma * rho_water * g)
    print(f"The excess pressure ΔP2 is determined by a gravity-capillary balance: ΔP2 = 2 * sqrt(σ * ρ_water * g)")
    print(f"ΔP2 = 2 * sqrt({sigma} * {rho_water} * {g}) = {delta_p2:.2f} Pa")

    # Jet speed from Bernoulli's principle
    v_jet2 = math.sqrt((2 * delta_p2) / rho_air)
    print(f"The jet speed V2 is calculated using Bernoulli's principle: V2 = sqrt(2 * ΔP2 / ρ_air)")
    print(f"V2 = sqrt(2 * {delta_p2:.2f} / {rho_air}) = {v_jet2:.2f} m/s\n")

    # --- Step 4: Final Answer ---
    print("--- Final Result ---")
    print(f"The calculated speed for the 2 mm bubble is approximately {round(v_jet1)} m/s.")
    print(f"The calculated speed for the 2 cm bubble is approximately {round(v_jet2)} m/s.")
    print(f"The corresponding speeds are ({round(v_jet1)}, {round(v_jet2)}).")

calculate_jet_speed()