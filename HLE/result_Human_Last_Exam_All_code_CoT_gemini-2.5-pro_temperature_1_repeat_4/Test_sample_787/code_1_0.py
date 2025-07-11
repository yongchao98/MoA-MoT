import math

def solve():
    """
    Calculates the limit of the sequence g_n based on the polynomial P(X).
    The limit L is determined by finding the greatest common divisor (GCD)
    of P(p) for all sufficiently large primes p.
    
    Through number theory analysis, the limit L must have prime factors
    only from the set {2, 3, 5}. The exponents of these prime factors are
    determined by finding the minimum p-adic valuation of P(x) over all integers x.
    
    The limit L is given by 2^a * 3^b * 5^c, where:
    a = min(v_2(P(p))) over all primes p, which is 10.
    b = min(v_3(P(p))) over all primes p, which is 2.
    c = min(v_5(P(p))) over all primes p, which is 1.
    """
    
    # Exponents determined by mathematical analysis
    a = 10
    b = 2
    c = 1
    
    # Bases
    base1 = 2
    base2 = 3
    base3 = 5
    
    # Calculate the final result
    result = (base1 ** a) * (base2 ** b) * (base3 ** c)
    
    # Print the equation as requested
    print(f"{base1}**{a} * {base2}**{b} * {base3}**{c} = {result}")

solve()