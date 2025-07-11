import collections

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    return dict(factors)

def main():
    """
    Main function to solve the problem.
    """
    a1, a2, a3, a4 = 8, 59, 79, 29

    print("Step 1: Calculate the three terms that constitute the sum S.")
    # Term 1: from p1=2, exponent a1=8
    term1 = a1 + 1
    print(f"Term 1 = (a1 + 1) = ({a1} + 1) = {term1}")

    # Term 2: from p2=29, exponent a2=59
    # Sum(1 + 2 + ... + (a2+1))
    term2 = (a2 + 1) * (a2 + 2) // 2
    print(f"Term 2 = (a2 + 1) * (a2 + 2) / 2 = ({a2} + 1) * ({a2} + 2) / 2 = {term2}")

    # Term 3: from p3=59, a3=79 and p4=79, a4=29
    # E_k is the sum for even indices, O_k is the sum for odd indices
    E3 = sum(a3 - b3 + 1 for b3 in range(0, a3 + 1, 2))
    O3 = sum(a3 - b3 + 1 for b3 in range(1, a3 + 1, 2))
    E4 = sum(a4 - b4 + 1 for b4 in range(0, a4 + 1, 2))
    O4 = sum(a4 - b4 + 1 for b4 in range(1, a4 + 1, 2))
    term3 = E3 * E4 + O3 * O4
    print(f"Term 3 (calculated from sums over exponents of 59 and 79) = {term3}")
    
    # Calculate the total sum S
    # S = term1 * term2 * term3

    print("\nStep 2: Find the prime factorization of S by combining factorizations of the terms.")
    factors1 = get_prime_factorization(term1)
    factors2 = get_prime_factorization(term2)
    factors3 = get_prime_factorization(term3)
    
    print(f"Prime factorization of Term 1 ({term1}): {factors1}")
    print(f"Prime factorization of Term 2 ({term2}): {factors2}")
    print(f"Prime factorization of Term 3 ({term3}): {factors3}")

    s_factors = collections.defaultdict(int)
    for p, e in factors1.items():
        s_factors[p] += e
    for p, e in factors2.items():
        s_factors[p] += e
    for p, e in factors3.items():
        s_factors[p] += e

    s_factors = dict(sorted(s_factors.items()))
    print(f"Combined prime factorization of S: {s_factors}")

    print("\nStep 3: Calculate the number of divisors of S from its prime factorization.")
    
    num_divisors = 1
    equation_parts = []
    for p, e in s_factors.items():
        num_divisors *= (e + 1)
        equation_parts.append(f"({e} + 1)")
    
    equation_str = " * ".join(equation_parts)
    print(f"The number of divisors is the product of (exponent + 1) for each prime factor.")
    print(f"Equation: {equation_str} = {num_divisors}")
    print(f"\nThe final answer is the number of divisors of the sum, which is {num_divisors}.")

if __name__ == "__main__":
    main()
