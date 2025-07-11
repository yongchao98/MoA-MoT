import math

def get_prime_factorization_exponents(n):
    """
    Calculates the exponents of the prime factorization of n.
    For n = p1^e1 * p2^e2 * ..., returns [e1, e2, ...].
    For n=1, returns an empty list.
    """
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
    if n == 1:
        return [], []

    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    
    primes = sorted(factors.keys())
    exponents = [factors[p] for p in primes]
    
    return primes, exponents

def solve_cardinality(l):
    """
    Calculates and prints the cardinalities |U_l| and |T_l|.
    """
    print(f"For l = {l}:")
    
    if l == 1:
        # Case l = 1
        cardinality = 1
        print("A) The cardinality |U_1| is 1.")
        print("B) For l=1, the prime factorization has no exponents e_i, so the expression for l > 1 does not apply.")
        return

    # Case l > 1
    primes, exponents = get_prime_factorization_exponents(l)
    
    # Calculate tau(l^2) = product of (2*e_i + 1)
    tau_l_squared = 1
    for e in exponents:
        tau_l_squared *= (2 * e + 1)
        
    cardinality = tau_l_squared - 1

    # Part A
    print(f"A) The cardinality |U_{l}| is {cardinality}.")
    
    # Part B
    # Build the expression string for l > 1
    factorization_str = " * ".join([f"{p}^{e}" for p, e in zip(primes, exponents)])
    expr_parts = [f"(2*{e}+1)" for e in exponents]
    expr_str = "*".join(expr_parts)
    
    print(f"B) For l = {l} = {factorization_str}, the exponents are {exponents}.")
    print(f"   The expression for |T_{l}| is {expr_str} - 1 = {cardinality}.")


# --- Example Usage ---
# You can change the value of `l` to solve for other integers.
l = 60
solve_cardinality(l)

print("\n" + "="*20 + "\n")

l = 1
solve_cardinality(l)

print("\n" + "="*20 + "\n")

# A case where l is a prime power
l = 9
solve_cardinality(l)
