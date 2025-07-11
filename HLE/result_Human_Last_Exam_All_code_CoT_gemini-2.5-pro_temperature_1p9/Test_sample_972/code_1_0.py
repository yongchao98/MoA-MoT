import math

def calculate_final_amplitude(initial_amplitude, slab_length, alpha):
    """
    Calculates the amplitude of an electromagnetic wave after passing through a slab
    with time-varying properties.

    Args:
        initial_amplitude (float): The initial amplitude of the wave (A).
        slab_length (float): The length of the slab in meters (L).
        alpha (float): The time-variation parameter in s^-1.
    """
    # Speed of light in vacuum (m/s)
    c0 = 299792458.0

    # The formula for the amplitude at x=L is A * exp(-alpha * L / c0)
    exponent = -alpha * slab_length / c0
    final_amplitude = initial_amplitude * math.exp(exponent)

    # Output the result and the equation with numbers substituted
    print("Final Amplitude Formula: A * exp(-alpha * L / c0)")
    print("Calculation:")
    print(f"{initial_amplitude} * exp(-{alpha} * {slab_length} / {c0}) = {final_amplitude}")
    
    return final_amplitude

if __name__ == '__main__':
    # Example parameters for the calculation
    # Initial amplitude of the electric field in V/m
    A = 10.0
    # Length of the slab in meters
    L = 0.5
    # Time-variation parameter alpha in s^-1
    alpha = 1.0e8

    print(f"Given Parameters:\nInitial Amplitude (A) = {A} V/m\nSlab Length (L) = {L} m\nAlpha = {alpha} 1/s\n")
    
    calculate_final_amplitude(initial_amplitude=A, slab_length=L, alpha=alpha)
