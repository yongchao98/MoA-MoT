import math

def get_prime_factorization(n):
    """Computes the prime factorization of n as a dictionary."""
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
       factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def multiply_factorizations(dict1, dict2):
    """Combines two dictionaries of prime factors."""
    res = dict1.copy()
    for p, e in dict2.items():
        res[p] = res.get(p, 0) + e
    return res

def main():
    """
    Solves the problem by calculating the number of divisors of S.
    S = sum_{d|N} f(d) for N = 2^8 * 29^59 * 59^79 * 79^29.
    We derived that S = F(2^8) * F(29^59) * F(59^79 * 79^29).
    """

    # Component 1: F(2^8)
    # F(2^8) = 8 + 1 = 9
    F2_val = 9
    F2_factors = get_prime_factorization(F2_val)

    # Component 2: F(29^59)
    # F(29^59) = sum_{j=1}^{59+1} j = 60 * 61 / 2 = 1830
    F29_val = (60 * 61) // 2
    F29_factors = get_prime_factorization(F29_val)

    # Component 3: F(59^79 * 79^29)
    # This equals A_e * B_e + A_o * B_o where A corresponds to 59^79 and B to 79^29
    
    # Sums for the prime 59 with exponent 79:
    # A_e = sum of (80-i) for i even in [0, 79] -> sum_{k=0}^{39} (80-2k)
    Ae = 40 * 80 - 2 * (39 * 40 // 2) # = 3200 - 1560 = 1640
    # A_o = sum of (80-i) for i odd in [0, 79] -> sum_{k=0}^{39} (80-(2k+1))
    Ao = 40 * 79 - 2 * (39 * 40 // 2) # = 3160 - 1560 = 1600
    
    # Sums for the prime 79 with exponent 29:
    # B_e = sum of (30-j) for j even in [0, 29] -> sum_{k=0}^{14} (30-2k)
    Be = 15 * 30 - 2 * (14 * 15 // 2) # = 450 - 210 = 240
    # B_o = sum of (30-j) for j odd in [0, 29] -> sum_{k=0}^{14} (30-(2k+1))
    Bo = 15 * 29 - 2 * (14 * 15 // 2) # = 435 - 210 = 225
    
    # We factor the sum F(59^79 * 79^29) = A_e*B_e + A_o*B_o by hand to avoid large numbers.
    # A_e*B_e + A_o*B_o = 1640*240 + 1600*225 = 393600 + 360000 = 753600.
    # Prime factorization: 753600 = 2^6 * 3^1 * 5^2 * 157^1
    F59_79_factors = {2: 6, 3: 1, 5: 2, 157: 1}

    # Total S is the product of the three components. Combine their prime factorizations.
    S_factors = multiply_factorizations(F2_factors, F29_factors)
    S_factors = multiply_factorizations(S_factors, F59_79_factors)

    # Calculate the number of divisors from the final exponents
    num_divisors = 1
    equation_parts = []
    sorted_primes = sorted(S_factors.keys())
    
    print("The sum S is a product of three main components:")
    print(f"F(2^8) = {F2_val}")
    print(f"F(29^59) = {F29_val}")
    print(f"F(59^79 * 79^29) = {Ae*Be + Ao*Bo}")
    print()

    print("The prime factorization of the sum S has the following exponents:")
    for p in sorted_primes:
        exponent = S_factors[p]
        num_divisors *= (exponent + 1)
        equation_parts.append(f"({exponent} + 1)")
        print(f"  Prime {p}: exponent {exponent}")

    print()
    print("The number of divisors is the product of (exponent + 1) for each prime factor:")
    final_equation = " * ".join(equation_parts)
    print(f"Number of divisors = {final_equation} = {num_divisors}")

if __name__ == "__main__":
    main()