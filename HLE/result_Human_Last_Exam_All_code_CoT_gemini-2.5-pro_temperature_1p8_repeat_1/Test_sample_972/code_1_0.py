import math

def calculate_final_amplitude(initial_amplitude, alpha, L):
    """
    Calculates the amplitude of an electromagnetic wave after passing through
    a slab with time-varying electromagnetic properties.

    The formula for the final amplitude A_L is A * exp(-alpha * L / (2*c)),
    where:
    A = initial amplitude of the wave
    alpha = rate of change of the material properties (from ε_r = μ_r = αt + β)
    L = length of the slab
    c = speed of light in vacuum

    Args:
        initial_amplitude (float): The amplitude of the wave at the entrance (A).
        alpha (float): The time-variation parameter alpha.
        L (float): The length of the slab.

    Returns:
        float: The amplitude of the wave at the rightmost boundary.
    """
    # Speed of light in meters per second
    c = 299792458.0

    # Calculate the exponent term
    exponent = -alpha * L / (2 * c)

    # Calculate the final amplitude
    final_amplitude = initial_amplitude * math.exp(exponent)
    
    # Print the equation with all the numbers
    print("Final Amplitude Calculation:")
    print(f"A_L = A * exp(-α * L / (2 * c))")
    print(f"A_L = {initial_amplitude} * exp(-{alpha} * {L} / (2 * {c}))")
    print(f"A_L = {final_amplitude}")
    
    return final_amplitude

if __name__ == '__main__':
    # Example usage with placeholder values.
    # The user should replace these with their specific values.
    
    # Initial amplitude of the incoming wave (e.g., in V/m)
    A_in = 1.0
    
    # Time-variation parameter alpha (in 1/s).
    # This must be small compared to the wave's angular frequency for the model to be valid.
    alpha_param = 1e9 # Corresponds to a fast-changing medium
    
    # Length of the slab (in meters)
    L_slab = 0.1 

    # Calculate and print the final amplitude
    calculate_final_amplitude(A_in, alpha_param, L_slab)