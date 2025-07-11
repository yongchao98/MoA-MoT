import math

def get_prime_factors(n):
    """
    Finds the prime factors of a given integer n.
    Assumes n is squarefree as per the problem description.
    """
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return list(factors)

def calculate_bounds(N):
    """
    Calculates the upper bound for the maximum norm (k_k,inf)
    in relation to the covolume (V) for a given squarefree N.
    """
    if N <= 0:
        print("Error: N must be a positive integer.")
        return

    # Step 1: Find the prime factors of N.
    prime_factors = get_prime_factors(N)
    
    # Check if N is squarefree by verifying prime factors.
    product_of_factors = 1
    for p in prime_factors:
        product_of_factors *= p
    if product_of_factors != N:
        print(f"Error: The number {N} is not squarefree.")
        return

    print(f"For the squarefree natural number N = {N} with prime factors {prime_factors}:")

    # Step 2: Calculate the product term shared by both formulas.
    prod_p_plus_1 = 1
    for p in prime_factors:
        prod_p_plus_1 *= (p + 1)

    # Step 3: Calculate the upper bound for the maximum norm, k_k,inf.
    # The formula is U = (4/3) * product_{p|N}(p+1)
    upper_bound_k = (4.0 / 3.0) * prod_p_plus_1

    # Step 4: Calculate the covolume V.
    # The formula is V = (pi/3) * product_{p|N}(p+1)
    covolume_V = (math.pi / 3.0) * prod_p_plus_1

    # Step 5: Establish the relationship k_k,inf <= (4/pi) * V and print the results.
    constant_factor = 4.0 / math.pi
    
    print("\nThe calculated upper bound for the maximum norm is:")
    print(f"k_k,inf <= {upper_bound_k:.4f}")

    print("\nThe calculated covolume is:")
    print(f"V = {covolume_V:.4f}")

    print("\nThe relationship between the bound and the covolume is k_k,inf <= (4/pi) * V.")
    print("Verifying this relationship with the calculated values:")
    
    # Print the final equation with all numbers.
    print(f"\nFinal Equation:")
    print(f"k_k,inf <= (4/pi) * V")
    print(f"{upper_bound_k:.4f} <= ({constant_factor:.4f}) * ({covolume_V:.4f})")
    
    # The result of the right-hand side calculation
    rhs_value = constant_factor * covolume_V
    print(f"{upper_bound_k:.4f} <= {rhs_value:.4f}")

# Example with a squarefree natural number N = 10 (2 * 5)
N = 10
calculate_bounds(N)
