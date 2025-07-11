import collections

def get_prime_factorization(n):
    """
    Returns a dictionary representing the prime factorization of n.
    Keys are prime factors, values are their exponents.
    """
    factors = collections.defaultdict(int)
    num = n
    d = 2
    while d * d <= num:
        while num % d == 0:
            factors[d] += 1
            num //= d
        d += 1
    if num > 1:
        factors[num] += 1
    return factors

def main():
    """
    Solves the problem by calculating the sum S and then finding the number of its divisors.
    """
    print("This script calculates the number of divisors of the sum S.")
    
    # Step 1: Calculate the sum over the exponent b.
    # Sum_{b=0 to 59} (60-b) is equivalent to Sum_{i=1 to 60} i.
    sum_b = sum(i for i in range(1, 61))

    # Step 2: Calculate the sums over the exponent c.
    # Exponent c ranges from 0 to 79.
    sum_c_even = sum(80 - c for c in range(0, 80, 2))
    sum_c_odd = sum(80 - c for c in range(1, 80, 2))

    # Step 3: Calculate the sums over the exponent d.
    # Exponent d ranges from 0 to 29.
    sum_d_even = sum(30 - d for d in range(0, 30, 2))
    sum_d_odd = sum(30 - d for d in range(1, 30, 2))

    # Step 4: Calculate the combined sum over c and d.
    sum_cd = sum_c_even * sum_d_even + sum_c_odd * sum_d_odd
    
    # Step 5: Calculate the total sum S.
    S = 9 * sum_b * sum_cd
    
    print("\n--- Intermediate Calculations ---")
    print(f"Sum over b's term (1 to 60): {sum_b}")
    print(f"Sum over c's term (c even): {sum_c_even}")
    print(f"Sum over c's term (c odd): {sum_c_odd}")
    print(f"Sum over d's term (d even): {sum_d_even}")
    print(f"Sum over d's term (d odd): {sum_d_odd}")
    print(f"Total sum over c and d terms: {sum_cd}")
    print(f"Final value of S = 9 * {sum_b} * {sum_cd} = {S}")

    # Step 6: Find the prime factorization of S by factoring its components.
    factors_S = collections.defaultdict(int)
    
    # Factor each component of S
    components = {9, sum_b, sum_cd}
    for comp in components:
        factors_comp = get_prime_factorization(comp)
        for p, e in factors_comp.items():
            factors_S[p] += e
            
    print("\n--- Prime Factorization of S ---")
    s_fact_str = " * ".join([f"{p}^{e}" for p, e in sorted(factors_S.items())])
    print(f"S = {s_fact_str}")

    # Step 7: Calculate the number of divisors of S.
    num_divisors = 1
    equation_parts = []
    # The number of divisors is the product of (exponent + 1) for each prime factor.
    for p in sorted(factors_S.keys()):
        exponent = factors_S[p]
        num_divisors *= (exponent + 1)
        equation_parts.append(f"({exponent} + 1)")
        
    equation = " * ".join(equation_parts)
    
    print("\n--- Final Answer ---")
    print("The number of divisors of S is calculated from the exponents of its prime factors:")
    print(f"Number of divisors = {equation} = {num_divisors}")

if __name__ == "__main__":
    main()