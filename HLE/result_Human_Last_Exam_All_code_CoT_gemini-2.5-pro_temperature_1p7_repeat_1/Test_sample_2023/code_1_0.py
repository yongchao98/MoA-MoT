import math

def solve():
    """
    Calculates the value of f(alpha_p, beta_p, gamma_p) mod p for p = 2^127 - 1.
    """
    # p is the Mersenne prime 2^127 - 1
    p = (2**127) - 1

    # As derived in the explanation, the problem simplifies to calculating
    # -V^2 mod p, where V is the product of three multinomial coefficients.

    # V0 = C(6; 1, 4, 1) = 6! / (1! * 4! * 1!)
    V0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))

    # V1 = C(8; 3, 2, 3) = 8! / (3! * 2! * 3!)
    V1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))

    # V2 = C(10; 4, 2, 4) = 10! / (4! * 2! * 4!)
    V2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

    # V is the product of these values
    V = V0 * V1 * V2

    # The final expression to compute is -V^2 mod p.
    # In modular arithmetic, -x is equivalent to p-x.
    V_squared_mod_p = pow(V, 2, p)
    result = (p - V_squared_mod_p) % p
    
    # Printing the intermediate values as requested
    print(f"The calculation is performed modulo the prime p:")
    print(f"p = {p}")
    print("-" * 20)
    print(f"The calculation depends on three intermediate multinomial coefficients:")
    print(f"V0 = C(6; 1, 4, 1) = {V0}")
    print(f"V1 = C(8; 3, 2, 3) = {V1}")
    print(f"V2 = C(10; 4, 2, 4) = {V2}")
    print("-" * 20)
    print(f"These are multiplied to form the base V:")
    print(f"V = V0 * V1 * V2 = {V}")
    print("-" * 20)
    print(f"The final equation derived from the theory is:")
    print(f"Result = -(V^2) mod p")
    print("-" * 20)
    print(f"The final result is:")
    print(result)

solve()