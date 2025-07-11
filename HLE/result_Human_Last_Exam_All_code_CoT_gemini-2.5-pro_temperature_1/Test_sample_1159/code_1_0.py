import math

def solve_float_bits():
    """
    Calculates the fewest random bits required to generate a uniform
    random positive floating-point number.
    """
    # Parameters from the IEEE 754 64-bit standard example
    m = 52  # mantissa bits
    e = 11  # exponent bits

    print("This script calculates the fewest random bits to generate a uniform random positive float.")
    print(f"Using parameters: m = {m}, e = {e}\n")

    # Step 1: Count the total number of positive representable floating-point values.
    # This includes normal numbers, subnormal numbers, and zero.
    #
    # - Normal numbers: The exponent E can range from 1 to 2^e - 2.
    #   The number of such exponents is (2**e - 2).
    #   For each normal exponent, there are 2**m possible mantissas.
    num_normal_exponents = (2**e - 2)
    num_mantissas = 2**m
    num_normal_values = num_normal_exponents * num_mantissas

    # - Subnormal numbers: The exponent E is 0.
    #   There are 2**m possible mantissas, representing zero and all subnormal values.
    num_subnormal_values = num_mantissas

    # - Total positive values = (num_normal_values) + (num_subnormal_values)
    #   This simplifies to: total = (2**e - 2 + 1) * 2**m = (2**e - 1) * 2**m
    total_positive_values = (2**e - 1) * num_mantissas

    print("--- Calculation Steps ---")
    print("1. Counting all positive representable floating-point numbers:")
    print(f"   - Total Values = (Number of Normal Exponents + 1) * (Number of Mantissas)")
    print(f"   - Total Values = ((2^{e} - 2) + 1) * 2^{m}")
    print(f"   - Total Values = (2^{e} - 1) * 2^{m}\n")


    # Step 2: Calculate the number of bits required to uniquely select one value.
    # Bits = ceil(log2(total_positive_values))
    #      = ceil(log2((2**e - 1) * 2**m))
    #      = ceil(m + log2(2**e - 1))
    # Since e-1 < log2(2**e - 1) < e for e > 1, ceil(log2(2**e - 1)) = e.
    bits_for_exponent_part = math.ceil(math.log2(2**e - 1))
    total_bits = m + bits_for_exponent_part

    print("2. Calculating the minimum number of random bits required:")
    print(f"   - Bits = ceil(log2(Total Values))")
    print(f"   - Bits = m + ceil(log2(2^{e} - 1))")
    print(f"   - For e = {e}, ceil(log2(2^{e} - 1)) evaluates to {bits_for_exponent_part}")
    print(f"   - Thus, the symbolic formula is: m + e\n")
    
    print("--- Final Answer ---")
    print("The final equation for the number of bits is:")
    print(f"m + e = {m} + {e} = {total_bits}")

solve_float_bits()