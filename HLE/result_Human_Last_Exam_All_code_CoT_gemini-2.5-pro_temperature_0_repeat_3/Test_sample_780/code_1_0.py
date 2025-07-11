def solve():
    """
    This function calculates the value of S(N) mod 23627 based on the derived formula.
    """
    p = 23627
    k = 203
    k_minus_1 = 202

    # We need to compute (k * (k-1)^-1) mod p
    # This is equivalent to (k * pow(k-1, p-2, p)) mod p
    
    # Calculate the modular inverse of k-1
    # pow(base, exp, mod) is efficient for modular exponentiation
    inv_k_minus_1 = pow(k_minus_1, p - 2, p)
    
    # Calculate the final result
    result = (k * inv_k_minus_1) % p
    
    # Print the equation with all the numbers
    print(f"The problem simplifies to calculating (k / (k-1)) mod p.")
    print(f"The equation is: result = ({k} * pow({k_minus_1}, -1, {p})) % {p}")
    print(f"The modular inverse of {k_minus_1} modulo {p} is: {inv_k_minus_1}")
    print(f"So, the calculation is: ({k} * {inv_k_minus_1}) % {p}")
    print(f"The final result is: {result}")

solve()