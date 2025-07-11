import math

def solve_floating_point_bits():
    """
    Calculates the fewest random bits required to generate a uniform random
    n-bit floating point number over the unit interval [0, 1].
    """
    m_sym = 'm'
    e_sym = 'e'
    B_sym = 'B'

    print("This script calculates the fewest random bits required to generate a uniform random floating-point number in [0, 1].")
    print("\nThe steps are as follows:")

    print("\n1. Count the number of representable floating-point numbers in [0, 1] (N).")
    print("   - For a number to be in [0, 1], the sign bit 's' must be 0.")
    print("   - The count includes zero, subnormal numbers, and normalized numbers.")
    print(f"   - The total number of such representations, N, is a function of mantissa bits ({m_sym}) and bias ({B_sym}).")
    print(f"   - The derived count is: N = {B_sym} * 2**{m_sym} + 1")

    print(f"\n2. Relate the bias B to the number of exponent bits {e_sym}.")
    print(f"   - A common and mathematically convenient choice for the bias is {B_sym} = 2**({e_sym}-1).")
    print(f"   - Substituting this into the count N gives:")
    print(f"   N = (2**({e_sym}-1)) * 2**{m_sym} + 1 = 2**({m_sym} + {e_sym} - 1) + 1")

    print("\n3. Calculate the minimum number of random bits required.")
    print("   - This is given by the ceiling of the base-2 logarithm of N: ceil(log2(N)).")
    print(f"   - Bits = ceil(log2(2**({m_sym} + {e_sym} - 1) + 1))")
    print(f"   - The value log2(2**({m_sym}+{e_sym}-1) + 1) is strictly greater than {m_sym}+{e_sym}-1 and strictly less than {m_sym}+{e_sym}.")
    print(f"   - Therefore, the ceiling is {m_sym} + {e_sym}.")

    print("\n--- Final Derivation ---")
    print(f"bits_needed = ceil(log2( N ))")
    print(f"            = ceil(log2( {B_sym} * 2**{m_sym} + 1 ))")
    print(f"            = ceil(log2( (2**({e_sym}-1)) * 2**{m_sym} + 1 ))")
    print(f"            = ceil(log2( 2**({m_sym} + {e_sym} - 1) + 1 ))")
    print(f"            = {m_sym} + {e_sym}")

solve_floating_point_bits()