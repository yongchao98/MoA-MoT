import collections

def prime_factorize(n):
    """
    Returns a dictionary of the prime factorization of n.
    """
    factors = collections.defaultdict(int)
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def main():
    """
    Main function to solve the problem.
    """
    # Step 1 & 5: Define exponents and calculate components
    # N_odd = 29^59 * 59^79 * 79^29
    # p1=29 (1 mod 4), a1=59
    # p2=59 (3 mod 4), a2=79
    # p3=79 (3 mod 4), a3=29
    a1, a2, a3 = 59, 79, 29

    # Calculate F_tau(N_odd) = F_tau(p1^a1) * F_tau(p2^a2) * F_tau(p3^a3)
    # F_tau(p^a) = (a+1)(a+2)/2
    F_tau_p1 = (a1 + 1) * (a1 + 2) // 2
    F_tau_p2 = (a2 + 1) * (a2 + 2) // 2
    F_tau_p3 = (a3 + 1) * (a3 + 2) // 2
    F_tau_N_odd = F_tau_p1 * F_tau_p2 * F_tau_p3

    # Calculate F_h(N_odd) = F_h(p1^a1) * F_h(p2^a2) * F_h(p3^a3)
    # F_h(p^a) for p=1 mod 4 is (a+1)(a+2)/2
    F_h_p1 = (a1 + 1) * (a1 + 2) // 2
    # F_h(p^a) for p=3 mod 4 is floor(a/2)+1
    F_h_p2 = a2 // 2 + 1
    F_h_p3 = a3 // 2 + 1
    F_h_N_odd = F_h_p1 * F_h_p2 * F_h_p3

    # Step 4 & 6: Calculate S' and S
    # S' = 1/2 * (F_tau(N_odd) + F_h(N_odd))
    S_prime_numerator = F_tau_N_odd + F_h_N_odd
    S_prime = S_prime_numerator // 2
    
    # S = 9 * S'
    S = 9 * S_prime

    # Step 7: Find the number of divisors of S by first finding its prime factorization.
    # We factor the components of S to handle large numbers.
    # S = 9 * ( (F_tau_p1*F_tau_p2*F_tau_p3) + (F_h_p1*F_h_p2*F_h_p3) ) / 2
    # S = 9/2 * (F_h_p1 * (F_tau_p2*F_tau_p3 + F_h_p2*F_h_p3))
    # S = 9 * (F_h_p1/2) * (F_tau_p2*F_tau_p3 + F_h_p2*F_h_p3)
    
    term1 = 9
    term2 = F_h_p1 // 2
    term3 = F_tau_p2 * F_tau_p3 + F_h_p2 * F_h_p3
    
    factors_S = collections.defaultdict(int)

    # Combine factorizations of each term
    for p, e in prime_factorize(term1).items():
        factors_S[p] += e
    for p, e in prime_factorize(term2).items():
        factors_S[p] += e
    for p, e in prime_factorize(term3).items():
        factors_S[p] += e
        
    # Calculate the number of divisors from the combined factorization
    num_divisors = 1
    equation_parts = []
    
    # Sort primes for consistent output
    sorted_primes = sorted(factors_S.keys())
    
    exponent_list = [str(factors_S[p]) for p in sorted_primes]
    print(f"The prime factorization of S has the form p1^e1 * p2^e2 * ...")
    print(f"The primes are: {', '.join(map(str, sorted_primes))}")
    print(f"The exponents are: {', '.join(exponent_list)}")
    
    for p in sorted_primes:
        exponent = factors_S[p]
        num_divisors *= (exponent + 1)
        equation_parts.append(f"({exponent} + 1)")

    # Final print statement as requested
    final_equation = " * ".join(equation_parts)
    print("\nThe number of divisors is the product of (exponent + 1):")
    print(f"Number of divisors = {final_equation} = {num_divisors}")

if __name__ == "__main__":
    main()