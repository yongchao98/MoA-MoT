import math

def calculate_jet_speeds():
    """
    Calculates the gas jet speed from a bursting bubble at an air-water interface
    for two different bubble diameters, considering the physical regime of each bubble.
    """
    # --- 1. Define Physical Constants ---
    # Surface tension of water (N/m)
    gamma = 0.072
    # Density of air (kg/m^3)
    rho_air = 1.225
    # Density of water (kg/m^3)
    rho_water = 1000
    # Acceleration due to gravity (m/s^2)
    g = 9.8

    # --- Bubble radii in meters ---
    r1_m = (2 / 1000) / 2  # 2 mm diameter -> 0.001 m radius
    r2_m = (2 / 100) / 2   # 2 cm diameter -> 0.01 m radius

    # --- 2. Calculate Capillary Length ---
    capillary_length = math.sqrt(gamma / (g * rho_water))
    
    print("This script calculates the speed of gas jets from bursting bubbles of different sizes.")
    print("-" * 60)
    print(f"Key physical parameters used:")
    print(f"  Surface tension of water (gamma): {gamma} N/m")
    print(f"  Density of air (rho_air): {rho_air} kg/m^3")
    print(f"  Gravitational acceleration (g): {g} m/s^2")
    print(f"  Calculated capillary length of water: {capillary_length*1000:.2f} mm")
    print("-" * 60)

    # --- 3. Analyze the Small Bubble (2 mm diameter) ---
    print("Case 1: Bubble Diameter = 2 mm (Radius = 1 mm)")
    print("Regime: Radius (1 mm) < Capillary Length (2.7 mm) -> Small bubble (spherical cap)")
    print("The driving pressure is the Laplace pressure: Delta_P = 2 * gamma / R")
    v1 = math.sqrt((4 * gamma) / (rho_air * r1_m))
    print("Final Equation: V = sqrt( (4 * gamma) / (rho_air * R) )")
    print(f"Calculation: V = sqrt( (4 * {gamma}) / ({rho_air} * {r1_m}) )")
    print(f"Resulting jet speed: {v1:.1f} m/s\n")
    

    # --- 4. Analyze the Large Bubble (2 cm diameter) ---
    print("Case 2: Bubble Diameter = 2 cm (Radius = 10 mm)")
    print("Regime: Radius (10 mm) > Capillary Length (2.7 mm) -> Large bubble (gravity-flattened cap)")
    print("Driving pressure is limited by gravity: Delta_P = 2 * sqrt(gamma * g * rho_water)")
    delta_p_large = 2 * math.sqrt(gamma * g * rho_water)
    v2 = math.sqrt((2 * delta_p_large) / rho_air)
    print("Final Equation: V = sqrt( (2 * 2 * sqrt(gamma * g * rho_water)) / rho_air )")
    print(f"Calculation: V = sqrt( (2 * 2 * sqrt({gamma} * {g} * {rho_water})) / {rho_air} )")
    print(f"Resulting jet speed: {v2:.1f} m/s")
    print("-" * 60)
    
    print(f"The calculated speeds for 2 mm and 2 cm diameter bubbles are approximately {round(v1)} m/s and {round(v2)} m/s respectively.")

if __name__ == "__main__":
    calculate_jet_speeds()