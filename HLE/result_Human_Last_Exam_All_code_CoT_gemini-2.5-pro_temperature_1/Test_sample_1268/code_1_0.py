import math

def calculate_field_properties(N):
    """
    This function calculates and displays properties for the imaginary quadratic field Q(sqrt(-N)).

    The user is asking for the relationship between the covolume (V) and an upper bound
    for an ideal's norm (which we interpret as the Minkowski bound, B), for a field
    defined by a squarefree natural number N. This script demonstrates this relationship
    for an imaginary quadratic field.
    """
    print(f"--- Analysis for the Imaginary Quadratic Field K = Q(sqrt(-{N})) ---")

    # A squarefree natural number N defines the field K = Q(sqrt(-N)).
    # The discriminant of K depends on N modulo 4. Let d = -N.
    d = -N
    if d % 4 == 1:
        discriminant = d
    else:
        discriminant = 4 * d
    
    abs_discriminant = abs(discriminant)
    print(f"The discriminant of the field is: {discriminant}")

    # The covolume (V) of the ring of integers of K is sqrt(|discriminant|) / 2.
    covolume_V = math.sqrt(abs_discriminant) / 2
    print(f"The covolume V is: {covolume_V:.5f}")

    # The Minkowski bound (B), which we interpret as k_{k,infty}, is an upper bound
    # on the norm of an ideal needed to represent any ideal class.
    # For an imaginary quadratic field, B = (2/pi) * sqrt(|discriminant|).
    minkowski_bound_B = (2 / math.pi) * math.sqrt(abs_discriminant)
    print(f"The Minkowski bound B is: {minkowski_bound_B:.5f}")

    # The relationship between B and V is B = (4/pi) * V.
    # We will now print the final equation with the calculated values to verify this.
    constant_factor = 4 / math.pi
    
    print("\nThe relationship between the bound B and covolume V is: B = (4/pi) * V")
    print("\nPlugging in the calculated values, we get the final equation:")
    
    # Final output as requested: "output each number in the final equation!"
    print(f"{minkowski_bound_B:.5f} = {constant_factor:.5f} * {covolume_V:.5f}")


# Let's use the squarefree natural number N=5 as an example.
# You can change this value to any other squarefree natural number.
N = 5
calculate_field_properties(N)