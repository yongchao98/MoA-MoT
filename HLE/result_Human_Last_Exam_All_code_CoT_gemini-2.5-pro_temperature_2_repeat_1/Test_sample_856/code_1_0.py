import math

def get_prime_factorization_of_factorial(n):
    """
    Calculates the prime factorization of n! using Legendre's formula.
    Returns a dictionary of {prime: exponent}.
    """
    # First, find all primes up to n using a simple sieve
    primes = []
    is_prime = [True] * (n + 1)
    for p in range(2, n + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, n + 1, p):
                is_prime[i] = False

    factors = {}
    for p in primes:
        exponent = 0
        p_power = p
        while p_power <= n:
            exponent += n // p_power
            p_power *= p
        factors[p] = exponent
    return factors

def euler_totient_from_factors(factors_dict):
    """
    Calculates Euler's totient function phi(N) given the prime
    factorization of N as a dictionary.
    phi(p^k) = p^(k-1) * (p-1)
    """
    result = 1
    for prime, exponent in factors_dict.items():
        phi_pk = (prime ** (exponent - 1)) * (prime - 1)
        result *= phi_pk
    return result

def solve_manifold_problem():
    """
    Solves the problem of counting 3-manifolds with |pi_1| = 10!.
    """
    n = 10
    group_order = math.factorial(n)

    print("Step 1: Identify possible fundamental group structures.")
    print(f"The order of the fundamental group is {n}! = {group_order}.")
    print("Based on the classification of finite groups that can act freely on the 3-sphere, there are two possibilities for a group of this order:")
    print(f"- A cyclic group of order {group_order}.")
    print(f"- A generalized quaternionic group of order {group_order} (since {group_order} is divisible by 4).\n")
    
    print("Step 2: Count the number of manifolds for each group structure.")

    # Case 1: Generalized Quaternionic Group
    num_quaternionic = 1
    print(f"For the non-cyclic, generalized quaternionic group, there is exactly {num_quaternionic} unique manifold.")

    # Case 2: Cyclic Group (Lens Spaces)
    print("For the cyclic group, the number of non-homeomorphic manifolds (lens spaces) is phi(N) / 2, where N is the group order.")
    
    # Calculate phi(10!)
    factors_10_fact = get_prime_factorization_of_factorial(n)
    
    factor_str_parts = [f"{p}^{e}" for p, e in factors_10_fact.items()]
    factor_str = " * ".join(factor_str_parts)
    print(f"First, we find the prime factorization of {n}! = {factor_str}.")

    phi_val = euler_totient_from_factors(factors_10_fact)
    print(f"Using this, we calculate Euler's totient function: phi({group_order}) = {phi_val}.")
    
    # Since 10! > 2, we divide by 2
    num_lens_spaces = phi_val // 2
    print(f"The number of lens spaces is {phi_val} / 2 = {num_lens_spaces}.\n")

    # Step 3: Sum the results
    total_manifolds = num_quaternionic + num_lens_spaces
    print("Step 3: Calculate the total number of manifolds by summing the counts from both cases.")
    print("\n--- Final Calculation ---")
    print(f"Total Manifolds = (Number from quaternionic case) + (Number from cyclic case)")
    print(f"Total Manifolds = {num_quaternionic} + {num_lens_spaces}")
    print(f"{total_manifolds} = {num_quaternionic} + {num_lens_spaces}")

if __name__ == '__main__':
    solve_manifold_problem()