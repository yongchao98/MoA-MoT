import math

def calculate_final_amplitude(A, L, alpha):
    """
    Calculates the amplitude of an electric field at the boundary of a time-varying slab.

    The formula for the output amplitude (A_out) is derived from Maxwell's equations
    for the specific medium described:
    A_out = A * exp(-(alpha * L) / (2 * c))

    Args:
        A (float): Amplitude of the incoming wave (e.g., in V/m).
        L (float): Length of the slab (in meters).
        alpha (float): Time-variation parameter of the slab's properties (in 1/s).
    """
    # Physical constant: speed of light in vacuum (m/s)
    c = 299792458

    # Calculate the terms for the final equation
    exponent_numerator = alpha * L
    exponent_denominator = 2 * c
    exponent = -exponent_numerator / exponent_denominator
    exp_factor = math.exp(exponent)
    A_out = A * exp_factor

    # Print the explanation and the final equation with values
    print("The final amplitude (A_out) is calculated using the formula:")
    print("A_out = A * exp(-(alpha * L) / (2 * c))\n")

    print(f"Given Parameters:")
    print(f"  Initial Amplitude (A) = {A}")
    print(f"  Slab Length (L) = {L}")
    print(f"  Time Variation (alpha) = {alpha}")
    print(f"  Speed of Light (c) = {c}\n")

    print("Final Equation with numbers:")
    # The prompt requires outputting each number in the final equation.
    # The following print statements show the step-by-step substitution.
    print(f"A_out = {A} * exp(-({alpha} * {L}) / (2 * {c}))")
    print(f"A_out = {A} * exp(-({exponent_numerator}) / ({exponent_denominator}))")
    print(f"A_out = {A} * exp({exponent})")
    print(f"A_out = {A} * {exp_factor}")
    print(f"A_out = {A_out}\n")

    print("--- Final Result ---")
    print(f"The amplitude of the electric field at x=L is: {A_out}")
    
    # Return value to be captured for the final answer format
    return A_out

# --- User-configurable values ---
# You can change these values to see the result for different scenarios.
initial_amplitude = 10.0  # V/m
slab_length = 1.0       # meters
alpha_parameter = 1.0e8 # 1/s
# --- End of user-configurable values ---

# Execute the calculation and store the result
final_amplitude_value = calculate_final_amplitude(initial_amplitude, slab_length, alpha_parameter)
# The final answer format is not required for a function that prints its own result,
# but we can prepare it just in case. The physics result is the formula itself.
# In this case, the question is implicitly asking for the formula A * exp(-(alpha * L) / (2 * c))
# Let's output the final numerical answer as requested by the format.
# print(f'<<<{final_amplitude_value}>>>') # This is for internal logic, the visible output is above.