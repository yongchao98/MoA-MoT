import math

def calculate_stable_pore_radius():
    """
    Calculates the stable radius of a pore during sintering where the
    sintering driving pressure is balanced by internal gas pressure.
    
    The equation being solved is: 2*gamma/r = P_gas
    where P_gas = nRT/V and V = (4/3)*pi*r^3
    
    This simplifies to: r = sqrt((3 * n * R * T) / (8 * pi * gamma))
    """
    # --- Assumptions and Constants ---
    
    # gamma: Surface energy of the ceramic (e.g., Alumina). Unit: Joules/meter^2
    gamma = 1.0
    
    # T: Isothermal sintering temperature in Kelvin. (e.g., 1600 C = 1873.15 K)
    T_celsius = 1600
    T_kelvin = T_celsius + 273.15
    
    # n: Moles of trapped gas in a single pore. This is a very small number.
    n_moles = 1e-15 # mol
    
    # R: Ideal gas constant. Unit: Joules / (mol * Kelvin)
    R = 8.314
    
    # --- Calculation ---
    numerator = 3 * n_moles * R * T_kelvin
    denominator = 8 * math.pi * gamma
    
    # Calculate the square of the radius first
    r_squared = numerator / denominator
    
    # The final radius in meters
    r_meters = math.sqrt(r_squared)
    
    # Convert to a more readable unit (micrometers)
    r_micrometers = r_meters * 1e6
    
    # --- Output Results ---
    print("--- Sintering Pore Stabilization Calculation ---")
    print(f"This script calculates the stable radius of a pore where shrinkage is halted by trapped gas.")
    print("\nInput parameters based on physical constants and assumptions:")
    print(f"Surface Energy (gamma): {gamma} J/m^2")
    print(f"Sintering Temperature (T): {T_celsius}°C ({T_kelvin:.2f} K)")
    print(f"Moles of Trapped Gas (n): {n_moles:.1e} mol")
    print(f"Ideal Gas Constant (R): {R} J/(mol*K)")
    
    print("\nDerived Equation: r = sqrt((3 * n * R * T) / (8 * pi * gamma))")
    
    # Printing each number in the final equation calculation
    print(f"\nFinal calculation is: sqrt((3 * {n_moles:.1e} * {R} * {T_kelvin:.2f}) / (8 * {math.pi:.4f} * {gamma}))")

    print("\n--- Result ---")
    print(f"The pore will stop shrinking at a calculated radius of approximately {r_micrometers:.2f} micrometers.")
    print("This demonstrates how a small amount of trapped gas can lead to residual porosity.")

calculate_stable_pore_radius()
<<<D>>>