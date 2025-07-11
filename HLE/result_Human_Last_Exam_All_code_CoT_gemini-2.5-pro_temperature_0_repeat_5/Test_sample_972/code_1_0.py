import math

def solve_em_amplitude():
    """
    This function calculates and prints the symbolic formula for the amplitude
    of an electric field at the boundary of a time-varying slab.
    """

    # Symbolic variables from the problem description
    A = "A"  # Initial amplitude of the wave
    alpha = "α" # Parameter for time variation of the slab
    L = "L"  # Length of the slab
    c = "c"  # Speed of light in vacuum

    # The final derived formula for the amplitude at x=L is:
    # Amplitude(L) = A * exp(-(alpha * L) / (2 * c))

    # We print the final equation step-by-step to show all its components,
    # including the numbers.
    print("The amplitude of the electric field at the rightmost boundary (x=L) is given by the formula:")
    print(f"Amplitude({L}) = {A} * exp( [ -1 * ({alpha} * {L}) ] / [ 2 * {c} ] )")
    print("\nWhere:")
    print(f"  {A}: Initial amplitude of the electric field.")
    print(f"  {L}: The length of the slab.")
    print(f"  {alpha}: The parameter governing the time-variation of the slab's properties.")
    print(f"  {c}: The speed of light in vacuum.")
    print(f"  exp(): The exponential function.")
    print("The numbers in the final equation are:")
    print(f"  -1: The coefficient in the numerator of the exponent.")
    print(f"  2: The coefficient in the denominator of the exponent.")

solve_em_amplitude()
