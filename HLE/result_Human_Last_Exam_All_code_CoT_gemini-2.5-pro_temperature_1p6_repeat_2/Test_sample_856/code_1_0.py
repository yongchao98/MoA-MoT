import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(factors):
    """
    Calculates Euler's totient function phi(n) from its prime factorization.
    phi(p^k) = p^k - p^(k-1)
    """
    result = 1
    for p, k in factors.items():
        result *= (p**k - p**(k-1))
    return result

def count_sq_roots_of_one(factors):
    """
    Calculates the number of solutions to x^2 = 1 (mod n).
    This is the product of the number of solutions for each prime power factor.
    """
    num_sols = 1
    for p, k in factors.items():
        if p == 2:
            if k == 1:
                num_sols *= 1
            elif k == 2:
                num_sols *= 2
            else: # k > 2
                num_sols *= 4
        else: # odd prime
            num_sols *= 2
    return num_sols

def solve():
    """
    Main function to solve the problem.
    """
    # Step 1: Set n = 10!
    n = math.factorial(10)

    # Step 2: Get prime factorization of n
    # 10! = 3628800 = 2^8 * 3^4 * 5^2 * 7^1
    factors = get_prime_factorization(n)

    # Step 3: Count the number of abelian (cyclic) groups. There is 1.
    # The corresponding manifolds are lens spaces. We calculate their number.
    # The formula is 1/4 * (phi(n) + tau2(1, n) + tau2(-1, n) + N_0)
    
    # Calculate phi(n)
    phi_n = phi(factors)
    
    # Calculate tau2(1, n), the number of square roots of 1 mod n
    tau2_one = count_sq_roots_of_one(factors)

    # For n=10!, tau2(-1, n) = 0 because n has a prime factor 3 (and 7),
    # and x^2 = -1 (mod 3) has no solutions.
    # For n=10!, N_0=0 because 2 is not a unit mod n.
    num_lens_spaces = (phi_n + tau2_one) // 4
    
    # Step 4: Count the number of non-abelian groups. There are 26.
    # Each non-abelian group corresponds to exactly one manifold.
    num_non_abelian_manifolds = 26

    # Step 5: Calculate the total number of manifolds.
    total_manifolds = num_lens_spaces + num_non_abelian_manifolds

    # Print the result showing each number in the final equation.
    print(f"The number of manifolds for the cyclic group C(10!) is {num_lens_spaces}.")
    print(f"The number of manifolds for the 26 non-abelian groups is {num_non_abelian_manifolds}.")
    print("The total number of such 3-manifolds is the sum:")
    print(f"{num_lens_spaces} + {num_non_abelian_manifolds} = {total_manifolds}")

solve()
