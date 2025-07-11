def solve():
    """
    Calculates the S(n) mod p based on the derived formula.
    """
    p = 23627
    K = 203

    # The final result is given by the expression K * (K-1)^-1 mod p.
    # Let's compute the components of this expression.

    # 1. Calculate K-1
    val_to_invert = K - 1

    # 2. Calculate the modular multiplicative inverse of (K-1) mod p
    # using Python's pow(base, exponent, modulus) function.
    # For exponent = -1, it computes the modular inverse.
    inverse = pow(val_to_invert, -1, p)

    # 3. Calculate the final result
    result = (K * inverse) % p

    # Print out the final equation with all the numbers
    print(f"The final value is determined by the equation: (K * (K-1)^-1) mod p")
    print(f"With K = {K} and p = {p}:")
    final_equation = f"({K} * pow({val_to_invert}, -1, {p}))"
    print(f"The expression is: {final_equation} mod {p}")
    print(f"The inverse term pow({val_to_invert}, -1, {p}) = {inverse}")
    print(f"So the calculation becomes: ({K} * {inverse}) mod {p} = {result}")

solve()