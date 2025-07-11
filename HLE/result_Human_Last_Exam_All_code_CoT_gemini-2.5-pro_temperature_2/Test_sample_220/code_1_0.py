import math

def calculate_jet_speed():
    """
    Calculates the gas jet speed from a bursting bubble for two different diameters.
    The model choice depends on the bubble size relative to the capillary length.
    """

    # Physical constants
    sigma = 0.072  # Surface tension of water (N/m)
    rho_water = 1000  # Density of water (kg/m^3)
    rho_air = 1.2   # Density of air (kg/m^3)
    g = 9.8        # Acceleration due to gravity (m/s^2)

    # Bubble diameters
    d1 = 0.002  # 2 mm in meters
    d2 = 0.02   # 2 cm in meters

    r1 = d1 / 2
    r2 = d2 / 2

    # --- Calculation for the small bubble (D = 2 mm) ---
    # For small bubbles (R < capillary length), pressure is from surface tension (Laplace pressure).
    # v = sqrt(2 * ΔP / ρ_air) = sqrt(2 * (2 * σ / R) / ρ_air) = sqrt(4 * σ / (R * ρ_air))
    
    delta_p1 = 2 * sigma / r1
    v1 = math.sqrt(2 * delta_p1 / rho_air)
    
    print("--- For bubble diameter of 2 mm ---")
    print("Model: Surface tension dominated (Young-Laplace pressure).")
    print(f"Excess pressure ΔP = 2 * {sigma} N/m / {r1} m = {delta_p1:.2f} Pa")
    print("Jet speed v = sqrt(2 * ΔP / ρ_air)")
    print(f"v = sqrt(2 * {delta_p1:.2f} Pa / {rho_air} kg/m^3) = {v1:.2f} m/s")
    print("-" * 35)

    # --- Calculation for the large bubble (D = 2 cm) ---
    # For large bubbles (R > capillary length), pressure is from gravity flattening the bubble.
    # ΔP ≈ 2 * sqrt(σ * ρ_water * g)
    # v = sqrt(2 * ΔP / ρ_air)

    delta_p2 = 2 * math.sqrt(sigma * rho_water * g)
    v2 = math.sqrt(2 * delta_p2 / rho_air)

    print("--- For bubble diameter of 2 cm ---")
    print("Model: Gravity dominated (flattened bubble).")
    print(f"Excess pressure ΔP = 2 * sqrt({sigma} * {rho_water} * {g}) = {delta_p2:.2f} Pa")
    print("Jet speed v = sqrt(2 * ΔP / ρ_air)")
    print(f"v = sqrt(2 * {delta_p2:.2f} Pa / {rho_air} kg/m^3) = {v2:.2f} m/s")
    print("-" * 35)

    # Final rounded results
    final_v1 = round(v1)
    final_v2 = round(v2)
    
    print(f"The calculated speeds are approximately {final_v1} m/s for a 2 mm bubble and {final_v2} m/s for a 2 cm bubble.")

calculate_jet_speed()
<<<E>>>