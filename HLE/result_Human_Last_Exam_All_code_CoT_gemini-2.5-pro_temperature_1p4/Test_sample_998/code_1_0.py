import math
from collections import Counter

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary (Counter).
    """
    factors = Counter()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return factors

def main():
    """
    Solves the problem step-by-step.
    """
    # Let N = 2^a * 29^a1 * 59^b1 * 79^b2
    a = 8
    a1 = 59  # Exponent of 29 (1 mod 4)
    b1 = 79  # Exponent of 59 (3 mod 4)
    b2 = 29  # Exponent of 79 (3 mod 4)

    print("Step 1: The sum S is given by the formula:")
    print("S = sum_{m|N, m=1(mod 4)} tau(N/m)")
    print("This can be broken down into parts based on the exponents of the prime factors of N.\n")

    # Part 1: Contribution from the prime factor 2
    # For a divisor m = 2^x... to be 1 (mod 4), x must be 0.
    # The number of divisors of N/m = 2^(8-x)... is (8-x+1)...
    # Since x=0, this factor is (8-0+1) = 9.
    factor_2 = a + 1
    print(f"Step 2: Calculate the contribution from each prime factor.")
    print(f"For prime 2, exponent is {a}. The condition m=1(mod 4) requires the exponent of 2 in m to be 0.")
    print(f"This contributes a factor of (a+1) = ({a}+1) = {factor_2} to the sum S.\n")

    # Part 2: Contribution from the prime factor 29 (p = 1 mod 4)
    # This sum is over y from 0 to 59 for tau(29^(59-y)), which is (59-y+1).
    # sum_{y=0}^{a1} (a1+1-y) = sum_{k=1}^{a1+1} k
    sum_y = (a1 + 1) * (a1 + 2) // 2
    print(f"For prime 29, exponent is {a1}. The sum over its exponents is:")
    print(f"Sum_y = sum(y=0..{a1}) ({a1}+1-y) = {sum_y}\n")

    # Part 3: Contribution from primes 59 and 79 (q = 3 mod 4)
    # The condition is that for m = ...59^z * 79^w, z+w must be even.
    # Sum_{z,w} = sum_{z=0..b1, w=0..b2, z+w is even} (b1+1-z)(b2+1-w)
    
    # Let's calculate the sums for z and w separately for even/odd exponents
    # For z from 0 to b1=79
    sum_z_even = sum(b1 + 1 - z for z in range(0, b1 + 1) if z % 2 == 0)
    sum_z_odd = sum(b1 + 1 - z for z in range(0, b1 + 1) if z % 2 != 0)
    # For w from 0 to b2=29
    sum_w_even = sum(b2 + 1 - w for w in range(0, b2 + 1) if w % 2 == 0)
    sum_w_odd = sum(b2 + 1 - w for w in range(0, b2 + 1) if w % 2 != 0)
    
    sum_z_w = sum_z_even * sum_w_even + sum_z_odd * sum_w_odd
    print(f"For primes 59 (exp {b1}) and 79 (exp {b2}), the sum over their exponents (z, w) where z+w is even is:")
    print(f"Sum_z_w = (sum_z_even * sum_w_even) + (sum_z_odd * sum_w_odd)")
    print(f"Sum_z_w = ({sum_z_even} * {sum_w_even}) + ({sum_z_odd} * {sum_w_odd}) = {sum_z_w}\n")

    # Final value of S
    S = factor_2 * sum_y * sum_z_w
    print(f"Step 3: Calculate the total sum S.")
    print(f"S = {factor_2} * {sum_y} * {sum_z_w} = {S}\n")
    
    # Step 4: Find the number of divisors of S
    print(f"Step 4: Find the prime factorization of S.")
    s_factors = get_prime_factorization(S)
    
    factor_list = sorted(s_factors.items())
    factor_str = " * ".join([f"{p}^{e}" for p, e in factor_list])
    print(f"The prime factorization of S is: {factor_str}\n")
    
    print(f"Step 5: Calculate the number of divisors of S using the exponents.")
    num_divisors = 1
    eq_parts = []
    exp_parts = []
    for p, e in factor_list:
        num_divisors *= (e + 1)
        eq_parts.append(f"({e}+1)")
        exp_parts.append(str(e))

    # The problem asks to output each number in the final equation
    print("The number of divisors is the product of (exponent + 1) for each prime factor.")
    print(f"The exponents are: {', '.join(exp_parts)}")
    print(f"Number of divisors = {' * '.join(eq_parts)} = {num_divisors}")
    
    print("\nFinal Answer:")
    print(num_divisors)

if __name__ == '__main__':
    main()