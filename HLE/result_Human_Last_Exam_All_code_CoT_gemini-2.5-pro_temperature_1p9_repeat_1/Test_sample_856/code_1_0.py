import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors and their exponents for a given number.
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve_manifold_problem():
    """
    Solves the problem of counting 3-manifolds with |pi_1| = 10!
    by analyzing the prime factors of the group's order.
    """
    n = 10
    target_order = math.factorial(n)

    print("Step 1: State the relevant topological theorem.")
    print("A closed orientable 3-manifold with a finite fundamental group G must be a spherical space form S³/G.")
    print("A necessary condition for a finite group G to act freely on S³ is that it has non-trivial Sylow p-subgroups for at most ONE odd prime p.\n")

    print(f"Step 2: Analyze the order of the fundamental group, which is {n}! = {target_order}.")
    print("We will find the prime factorization of 10! to check this condition.\n")

    # Get the prime factorization of 10!
    factors = get_prime_factorization(target_order)

    # Prepare the factorization equation string and find odd primes
    equation_parts = []
    odd_prime_factors = []
    # Sort by prime number for a clean equation
    for prime in sorted(factors.keys()):
        exponent = factors[prime]
        equation_parts.append(f"{prime}^{exponent}")
        # A prime is odd if it's not 2
        if prime != 2:
            odd_prime_factors.append(prime)

    factorization_string = " * ".join(equation_parts)
    print("Step 3: The prime factorization of 10! is calculated.")
    print(f"The final equation for the factorization is: 10! = {factorization_string}\n")

    print("Step 4: Check the condition on the odd prime factors.")
    num_odd_primes = len(odd_prime_factors)
    print(f"The distinct odd prime factors in the factorization are: {odd_prime_factors}.")
    print(f"The number of distinct odd prime factors is {num_odd_primes}.\n")

    print("Step 5: Conclusion based on the condition.")
    print(f"Because the number of distinct odd prime factors ({num_odd_primes}) is greater than one, no group of order 10! can act freely on S³.")
    print("Therefore, no closed orientable 3-manifold can have a fundamental group of this order.")
    final_answer = 0
    print(f"\nThe number of such manifolds is {final_answer}.")


# Run the solver
solve_manifold_problem()
<<<0>>>