from sympy import bernoulli, Rational, pi, factorint

def get_zeta_rational_part(n):
    """
    Calculates the rational part C_n where zeta(2n) = C_n * pi^(2n).
    C_n = abs(B_2n) * 2^(2n-1) / (2n)!
    """
    if n == 0:
        return Rational(-1, 2)
    
    # bernoulli(2*n) is B_{2n}
    B_2n = bernoulli(2*n)
    
    numerator = abs(B_2n) * (2**(2*n - 1))
    
    denominator = 1
    for i in range(1, 2*n + 1):
        denominator *= i
        
    return numerator / denominator

def main():
    """
    Calculates the sum as zeta(6)*zeta(8)/zeta(24) and expresses it
    as a rational number times an integer power of pi.
    """
    # Our sum is zeta(6) * zeta(8) / zeta(24)
    
    # zeta(6) = C_6 * pi^6
    C6 = get_zeta_rational_part(3)
    
    # zeta(8) = C_8 * pi^8
    C8 = get_zeta_rational_part(4)
    
    # zeta(24) = C_24 * pi^24
    C24 = get_zeta_rational_part(12)
    
    # The final rational coefficient is C6 * C8 / C24
    final_rational_coeff = (C6 * C8) / C24
    
    # The final power of pi is 6 + 8 - 24 = -10
    pi_power = -10
    
    # Output the result
    print("The sum is given by the expression zeta(6) * zeta(8) / zeta(24).")
    print("zeta(6) = (pi^6) / 945")
    print("zeta(8) = (pi^8) / 9450")
    print("The final expression is a rational number times an integer power of pi.")
    print("\nFinal calculation:")
    # Using 'print' for each number in the final equation as requested
    print(f"({final_rational_coeff.p} / {final_rational_coeff.q}) * pi^({pi_power})")
    
main()
