import math
import scipy.constants

def calculate_output_amplitude(A, L, alpha):
    """
    Calculates the amplitude of an electromagnetic wave after passing through
    a slab with time-varying properties.

    The final amplitude E_out is given by the formula:
    E_out = A * exp(-(alpha * L) / (2 * c_0))
    where c_0 is the speed of light in vacuum.

    Args:
        A (float): Initial amplitude of the electric field (in V/m).
        L (float): Length of the slab (in meters).
        alpha (float): Rate of change of the medium's properties (in 1/s).
    """
    # Use the speed of light from scipy constants
    c_0 = scipy.constants.c

    # Calculate the exponent term
    exponent = -(alpha * L) / (2 * c_0)

    # Calculate the final amplitude
    E_out = A * math.exp(exponent)
    
    # Print the equation and the values used
    print("The final amplitude (E_out) is calculated using the formula:")
    print("E_out = A * exp(-(alpha * L) / (2 * c_0))")
    print("\nUsing the following values:")
    print(f"Initial Amplitude (A)       = {A:.4f} V/m")
    print(f"Slab Length (L)             = {L:.4f} m")
    print(f"Rate of Change (alpha)      = {alpha:.4f} 1/s")
    print(f"Speed of Light (c_0)        = {c_0:.4f} m/s")

    # Print the final numerical equation
    print("\nPlugging in the numbers:")
    print(f"E_out = {A:.4f} * exp(-({alpha:.4f} * {L:.4f}) / (2 * {c_0:.4f}))")
    
    # Print the final result
    print(f"\nThe calculated output amplitude is: {E_out:.8f} V/m")

if __name__ == '__main__':
    # Example values for the parameters
    initial_amplitude = 1.0  # V/m
    slab_length = 1000.0     # meters
    alpha_rate = 5.0e4     # 1/s
    
    calculate_output_amplitude(initial_amplitude, slab_length, alpha_rate)
