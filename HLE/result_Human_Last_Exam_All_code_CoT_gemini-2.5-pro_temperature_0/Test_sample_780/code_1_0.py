def solve():
    """
    Calculates S(N) mod p based on the derived formula S(0) = K * (K-1)^-1 mod p.
    """
    p = 23627
    K = 203
    
    # The value we need to find is S(N) mod p, which simplifies to S(0) mod p.
    # The formula for S(0) mod p is K * (K-1)^(-1) mod p.
    
    K_minus_1 = K - 1
    
    # We need to find the modular multiplicative inverse of (K-1) modulo p.
    # pow(a, -1, m) calculates (a^-1) mod m.
    try:
        inv_K_minus_1 = pow(K_minus_1, -1, p)
    except ValueError:
        print(f"The inverse of {K_minus_1} modulo {p} does not exist.")
        return

    # Calculate the final result
    result = (K * inv_K_minus_1) % p
    
    # Output the numbers used in the final equation
    print(f"p = {p}")
    print(f"K = {K}")
    print(f"The expression to calculate is: (K * (K-1)^-1) mod p")
    print(f"This becomes: ({K} * {K_minus_1}^-1) mod {p}")
    print(f"The modular inverse of {K_minus_1} mod {p} is: {inv_K_minus_1}")
    print(f"The final calculation is: ({K} * {inv_K_minus_1}) mod {p}")
    print(f"Result: {result}")

solve()