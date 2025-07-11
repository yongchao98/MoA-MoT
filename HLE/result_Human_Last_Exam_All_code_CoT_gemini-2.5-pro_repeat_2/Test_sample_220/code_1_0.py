import math

def calculate_jet_speed(diameter_m, sigma, rho_gas):
    """
    Calculates the gas jet speed from a bursting bubble based on the Laplace pressure model.

    Args:
        diameter_m (float): The diameter of the bubble in meters.
        sigma (float): The surface tension of the liquid in N/m.
        rho_gas (float): The density of the gas in kg/m^3.

    Returns:
        float: The calculated speed of the gas jet in m/s.
    """
    radius = diameter_m / 2.0
    
    # The model is v_jet = sqrt(4 * sigma / (rho_gas * R))
    # This comes from equating Laplace pressure (2*sigma/R) with jet dynamic pressure (0.5*rho_gas*v^2)
    
    term_inside_sqrt = (4 * sigma) / (rho_gas * radius)
    speed = math.sqrt(term_inside_sqrt)
    
    print(f"--- Calculation for {diameter_m*1000:.0f} mm bubble ---")
    print(f"Equation: v = sqrt(4 * σ / (ρ_gas * R))")
    print(f"v = sqrt(4 * {sigma} / ({rho_gas} * {radius}))")
    print(f"v = sqrt({4 * sigma} / {rho_gas * radius})")
    print(f"v = sqrt({term_inside_sqrt})")
    print(f"Calculated speed: {speed:.1f} m/s\n")
    
    return speed

# Define physical constants
sigma_water = 0.072  # Surface tension of water (N/m)
rho_air = 1.2      # Density of air (kg/m^3)

# Define bubble diameters
diameter1_m = 0.002  # 2 mm
diameter2_m = 0.020  # 2 cm

# Perform and print calculations
v1 = calculate_jet_speed(diameter1_m, sigma_water, rho_air)
v2 = calculate_jet_speed(diameter2_m, sigma_water, rho_air)
