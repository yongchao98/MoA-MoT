import numpy as np
import scipy.constants

def calculate_inductance_change(h, R1):
    """
    Calculates the change in mutual inductance per unit length for two circuits
    inside an ideal cylindrical magnetic concentrator.

    Args:
        h (float): Separation distance between wires in each circuit (in meters).
        R1 (float): Inner radius of the concentrator shell (in meters).
    """
    # Physical constant: permeability of free space (mu_0) in H/m
    mu_0 = scipy.constants.mu_0

    # The formula for the change in mutual inductance per unit length is:
    # delta_M_per_length = (mu_0 * h^2) / (2 * pi * R1^2)
    delta_M_per_length = (mu_0 * h**2) / (2 * np.pi * R1**2)

    # --- Output Results ---
    print("The theoretical expression for the change in mutual inductance per unit length (ΔM/L) is:")
    print("ΔM/L = (μ₀ * h²) / (2 * π * R₁²)")
    print("\nThis script will now calculate the result using example values.")
    
    print("\n--- Equation and Values ---")
    # The prompt requires printing each number in the final equation.
    # We construct a string that shows the full calculation with substituted values.
    # Note: h**2 and R1**2 are calculated inline for clarity.
    
    calc_string = f"ΔM/L = ({mu_0:.4e} H/m * ({h} m)²) / (2 * {np.pi:.5f} * ({R1} m)²)"
    print(calc_string)
    
    # Breaking it down further to show each number explicitly as requested
    print("\nWhere the values are:")
    print(f"μ₀ = {mu_0}")
    print(f"h = {h}")
    print(f"2 = {2}")
    print(f"π = {np.pi}")
    print(f"R₁ = {R1}")

    print("\n--- Final Result ---")
    print(f"The calculated change in mutual inductance per unit length is: {delta_M_per_length:.4e} H/m")


if __name__ == '__main__':
    # Example values for the calculation.
    # These can be changed to see how the result varies.
    # We assume d >> h and h << R1.
    h_circuit = 0.01  # 1 cm separation between wires in a circuit
    R1_concentrator = 0.1  # 10 cm inner radius of the concentrator

    calculate_inductance_change(h=h_circuit, R1=R1_concentrator)