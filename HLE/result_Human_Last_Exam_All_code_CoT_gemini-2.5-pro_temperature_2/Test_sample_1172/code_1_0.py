import math

def calculate_inductance_change(d, h, R1, R2, L):
    """
    Calculates the change in mutual inductance between two circuits
    when a cylindrical magnetic concentrator is added, in the limit d >> h.

    Args:
        d (float): Distance between the two circuits (m).
        h (float): Separation of wires within each circuit (m).
        R1 (float): Inner radius of the concentrator shell (m).
        R2 (float): Outer radius of the concentrator shell (m).
        L (float): Length of the circuits (m).
    """
    # Check if the assumptions are met for the approximation
    if not d > 10 * h: # Heuristic check for d >> h
        print("Warning: The approximation d >> h may not be valid for the given inputs.")
    if not R2 > R1:
        print("Error: Outer radius R2 must be greater than inner radius R1.")
        return

    # Constants
    mu_0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)

    # Calculate components of the formula
    mu_factor = mu_0 / (2 * math.pi)
    concentrator_factor = (R2 / R1)**2 - 1
    geometry_factor = (h / d)**2

    # Calculate the total change in mutual inductance (ΔM)
    delta_M = mu_factor * concentrator_factor * geometry_factor * L

    # Output the result, showing each number in the final equation
    print("Derived formula for the change in mutual inductance:")
    print("ΔM = (μ₀ / (2π)) * ((R₂/R₁)² - 1) * (h/d)² * L\n")

    print("Calculation with provided values:")
    # Print the equation with substituted values
    print(f"ΔM = ({mu_factor:.2e}) * (({R2}/{R1})² - 1) * ({h}/{d})² * {L}")
    
    # Show intermediate calculation steps for clarity
    print(f"ΔM = ({mu_factor:.2e}) * ({concentrator_factor:.4f}) * ({geometry_factor:.4e}) * {L}")

    # Print the final numerical result in Henrys
    print(f"\nFinal Result: ΔM = {delta_M:.4e} H")


# --- Example Usage ---
# You can change these values to match a specific problem.
# Geometric parameters for the circuits and concentrator (in meters)
distance_between_circuits = 1.0   # d
wire_separation_in_circuit = 0.05 # h
concentrator_inner_radius = 0.2   # R1
concentrator_outer_radius = 0.4   # R2
circuit_length = 10.0             # L

calculate_inductance_change(distance_between_circuits, wire_separation_in_circuit,
                            concentrator_inner_radius, concentrator_outer_radius,
                            circuit_length)
