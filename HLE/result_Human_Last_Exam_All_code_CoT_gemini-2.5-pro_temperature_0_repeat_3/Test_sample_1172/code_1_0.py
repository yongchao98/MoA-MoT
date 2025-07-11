import numpy as np

def calculate_inductance_change(h, d, R1, mu_0):
    """
    Calculates the change in mutual inductance per unit length between two circuits
    when a cylindrical magnetic concentrator is added.

    The calculation is based on the formula derived using the method of images for
    the limit where d >> h.
    ΔM = (8 * μ₀ * h² * R₁²) / (π * (d² + 4 * R₁²)²)

    Args:
        h (float): Separation of wires in each circuit (m).
        d (float): Separation between the centers of the two circuits (m).
        R1 (float): Inner radius of the concentrator shell (m).
        mu_0 (float): Permeability of free space (H/m).

    Returns:
        float: The change in mutual inductance per unit length (H/m).
    """
    # Check if the geometric conditions for the approximation are met.
    if not d > 5 * h:  # A more concrete check for d >> h
        print(f"Warning: The approximation d >> h may not be very accurate as d/h = {d/h:.2f} is not very large.")
    if not d / 2 < R1:
        print(f"Warning: The circuits (at d/2={d/2} m) are not fully inside the inner radius of the shell (R1={R1} m).")

    # Print the symbolic formula
    print("The expression for the change in mutual inductance per unit length (ΔM) is:")
    print("ΔM = (8 * μ₀ * h² * R₁²) / (π * (d² + 4 * R₁²)²)")
    print("\nSubstituting the given values into the equation:")
    print(f"ΔM = (8 * {mu_0:.4e} H/m * ({h} m)² * ({R1} m)²) / (π * (({d} m)² + 4 * ({R1} m)²)²)")

    # Perform the calculation
    numerator = 8 * mu_0 * h**2 * R1**2
    denominator = np.pi * (d**2 + 4 * R1**2)**2
    delta_M = numerator / denominator
    
    return delta_M

if __name__ == '__main__':
    # Physical constant
    mu_0_val = 4 * np.pi * 1e-7  # Permeability of free space in H/m

    # --- User-defined Parameters ---
    # h: separation of wires in each circuit (meters)
    h_val = 0.01
    # d: separation between the centers of the two circuits (meters)
    d_val = 0.1
    # R1: inner radius of the concentrator shell (meters)
    R1_val = 0.2
    # -----------------------------

    # Calculate and print the result
    delta_M_result = calculate_inductance_change(h_val, d_val, R1_val, mu_0_val)
    
    print("\nCalculated change in mutual inductance per unit length:")
    print(f"ΔM = {delta_M_result:.4e} H/m")
