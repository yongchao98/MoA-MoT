import math

def calculate_jet_speeds():
    """
    Calculates the speed of gas jets from a bursting bubble at an air-water interface
    for two different bubble sizes, considering the dominant physical forces for each size.
    """
    # Physical constants
    sigma = 0.072  # Surface tension of water (N/m)
    rho_water = 1000  # Density of water (kg/m^3)
    rho_air = 1.2  # Density of air (kg/m^3)
    g = 9.8  # Acceleration due to gravity (m/s^2)

    # Bubble diameters
    d1_mm = 2
    d2_cm = 2
    
    # Convert diameters to radii in meters
    r1 = (d1_mm / 1000) / 2  # 2 mm diameter -> 0.001 m radius
    r2 = (d2_cm / 100) / 2  # 2 cm diameter -> 0.01 m radius

    # --- Calculation for the 2 mm bubble (surface tension dominates) ---
    print("--- For the 2 mm bubble (surface tension dominated regime) ---")
    
    # Young-Laplace equation for excess pressure in a hemispherical bubble
    delta_p1 = (2 * sigma) / r1
    
    # Bernoulli's principle to find jet velocity
    v1 = math.sqrt((2 * delta_p1) / rho_air)
    
    print(f"The radius is {r1} m.")
    print("The excess pressure is calculated using the Young-Laplace equation: ΔP = 2 * σ / R")
    print(f"ΔP = 2 * {sigma} / {r1} = {delta_p1:.2f} Pa")
    print("The jet velocity is calculated using Bernoulli's principle: v = sqrt(2 * ΔP / ρ_air)")
    print(f"v = sqrt(2 * {delta_p1:.2f} / {rho_air}) = {v1:.2f} m/s")
    print("-" * 50)

    # --- Calculation for the 2 cm bubble (gravity dominates) ---
    print("--- For the 2 cm bubble (gravity dominated regime) ---")
    
    # Capillary length, where gravity and surface tension effects are comparable
    capillary_length = math.sqrt(sigma / (rho_water * g))
    
    # For large bubbles, pressure is hydrostatic, related to the height (approx. 2 * capillary_length)
    # ΔP = ρ_water * g * H = ρ_water * g * (2 * l_c) = 2 * sqrt(σ * ρ_water * g)
    delta_p2 = 2 * math.sqrt(sigma * rho_water * g)
    
    # Bernoulli's principle to find jet velocity
    v2 = math.sqrt((2 * delta_p2) / rho_air)
    
    print(f"The radius is {r2} m.")
    print(f"The capillary length is l_c = sqrt(σ / (ρ_water * g)) = {capillary_length:.4f} m.")
    print("The bubble is large, so pressure is hydrostatic: ΔP = 2 * sqrt(σ * ρ_water * g)")
    print(f"ΔP = 2 * sqrt({sigma} * {rho_water} * {g}) = {delta_p2:.2f} Pa")
    print("The jet velocity is calculated using Bernoulli's principle: v = sqrt(2 * ΔP / ρ_air)")
    print(f"v = sqrt(2 * {delta_p2:.2f} / {rho_air}) = {v2:.2f} m/s")
    print("-" * 50)
    
    print(f"The calculated speeds are approximately {round(v1)} m/s for the 2 mm bubble and {round(v2)} m/s for the 2 cm bubble.")

calculate_jet_speeds()