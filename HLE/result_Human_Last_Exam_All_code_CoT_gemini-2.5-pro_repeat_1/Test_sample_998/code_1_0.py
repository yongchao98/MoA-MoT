import math
from collections import Counter

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary (Counter).
    """
    factors = Counter()
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def main():
    """
    This script solves the problem by calculating the value S and then finding the number of its divisors.
    Let N = 2^8 * 29^59 * 59^79 * 79^29.
    S = sum_{k|N} f(k), where f(k) is the number of divisors of k of the form 4m+1.
    This simplifies to S = sum_{d|N, d=1(mod 4)} d(N/d).
    A divisor d of N is 1 (mod 4) if d = 29^b1 * 59^b2 * 79^b3 with b2+b3 even.
    The term d(N/d) = 9 * (60-b1) * (80-b2) * (30-b3).
    S can be written as 9 * Sum1 * Sum2.
    """
    
    # Exponents from N
    exp_29 = 59
    exp_59 = 79
    exp_79 = 29

    # Calculate Sum1 = sum_{b1=0 to 59} (60-b1)
    # This is the sum of integers from 1 to 60.
    sum1 = (exp_29 + 1) * (exp_29 + 2) // 2

    # Calculate Sum2, which is sum_{b2,b3} (80-b2)(30-b3) where b2+b3 is even.
    # This splits into two cases: b2, b3 both even, or b2, b3 both odd.
    # Sum2 = (sum_{b2 even}(80-b2))*(sum_{b3 even}(30-b3)) + (sum_{b2 odd}(80-b2))*(sum_{b3 odd}(30-b3))

    # A = sum_{b2=0, even}^{79} (80-b2)
    A = sum(exp_59 + 1 - b2 for b2 in range(exp_59 + 1) if b2 % 2 == 0)
    
    # B = sum_{b3=0, even}^{29} (30-b3)
    B = sum(exp_79 + 1 - b3 for b3 in range(exp_79 + 1) if b3 % 2 == 0)

    # C = sum_{b2=0, odd}^{79} (80-b2)
    C = sum(exp_59 + 1 - b2 for b2 in range(exp_59 + 1) if b2 % 2 != 0)

    # D = sum_{b3=0, odd}^{29} (30-b3)
    D = sum(exp_79 + 1 - b3 for b3 in range(exp_79 + 1) if b3 % 2 != 0)
    
    sum2 = A * B + C * D
    
    # Now we find the prime factorization of S = 9 * sum1 * sum2
    s_factors = Counter({3: 2}) # from the factor 9
    
    sum1_factors = get_prime_factorization(sum1)
    s_factors.update(sum1_factors)
    
    sum2_factors = get_prime_factorization(sum2)
    s_factors.update(sum2_factors)

    # The exponents of the prime factors of S
    exponents = sorted(s_factors.values(), reverse=True)
    
    # Calculate the number of divisors
    num_divisors = 1
    for exp in exponents:
        num_divisors *= (exp + 1)
        
    # Format the final equation string
    equation_parts = [f"({exp}+1)" for exp in exponents]
    equation_str = " * ".join(equation_parts)
    
    print(f"The number of divisors is found from the exponents of the prime factorization of the sum S.")
    print(f"The exponents are {', '.join(map(str, exponents))}.")
    print(f"The number of divisors is {equation_str} = {num_divisors}")

if __name__ == "__main__":
    main()