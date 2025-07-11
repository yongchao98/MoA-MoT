import numpy as np

def calculate_final_amplitude(A, alpha, L):
    """
    Calculates the final amplitude of an EM wave passing through a time-varying slab.

    Args:
        A (float): Initial amplitude of the electric field (e.g., in V/m).
        alpha (float): Time-variation parameter of the slab (in 1/s).
        L (float): Length of the slab (in m).
    """
    # The speed of light in vacuum (m/s)
    c = 299792458.0

    # The derived formula for the final amplitude is A_final = A * exp(-alpha * L / (2 * c))

    # Calculate the value of the exponent
    exponent_val = -alpha * L / (2 * c)

    # Calculate the final amplitude
    final_amplitude = A * np.exp(exponent_val)

    # --- Output the results ---
    print("The formula for the final amplitude (A_final) is:")
    print("A_final = A * exp(-alpha * L / (2 * c))")
    print("\n--- Calculation with provided values ---")
    print(f"Initial Amplitude (A)      = {A} V/m")
    print(f"Slab parameter (alpha)     = {alpha:.2e} 1/s")
    print(f"Slab length (L)            = {L} m")
    print(f"Speed of light (c)         = {c:.2e} m/s")
    
    print("\nBreaking down the calculation of the exponent:")
    print(f"alpha * L                  = {alpha * L:.4e}")
    print(f"2 * c                      = {2 * c:.4e}")
    print(f"exponent = -(alpha*L)/(2*c) = {exponent_val:.4f}")

    print("\n--- Final Result ---")
    print(f"The final amplitude at x=L is: {final_amplitude:.6f} V/m")


if __name__ == '__main__':
    # Example usage with sample values
    initial_amplitude = 1.0  # V/m
    slab_alpha = 1.5e10      # 1/s
    slab_length = 0.1        # m (10 cm)
    
    calculate_final_amplitude(initial_amplitude, slab_alpha, slab_length)