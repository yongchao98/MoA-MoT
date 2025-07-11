import math
from collections import Counter

def get_prime_factorization(num):
    """Computes the prime factorization of a given number."""
    num = abs(int(num))
    if num <= 1:
        return Counter()
    
    factors = Counter()
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def solve_and_print():
    """
    This function solves the problem by calculating the eigenvalues of the Seidel matrix
    and then finding the LCM of their absolute values, which corresponds to the
    maximum order in the Smith group.
    """
    
    # Step 1: Define graph parameters
    n_set = 50
    k_subset = 5
    # Adjacency is for intersection size 3. For the eigenvalue formula, this
    # corresponds to a parameter t_dist = k_subset - 3 = 2.
    t_dist = k_subset - 3

    # Step 2: Calculate eigenvalues of the adjacency matrix A (theta_j)
    theta = {}
    for j in range(k_subset + 1):
        val = 0
        for i in range(t_dist + 1):
            term = ((-1)**i) * math.comb(j, i) * math.comb(k_subset - j, t_dist - i) * math.comb(n_set - k_subset - j, t_dist - i)
            val += term
        theta[j] = val
        
    # Step 3: Calculate eigenvalues of the Seidel matrix S (s_j)
    N = math.comb(n_set, k_subset)
    k_degree = theta[0]
    s = {}
    # Eigenvalue for the all-ones vector
    s[0] = N - 1 - 2 * k_degree
    # Other eigenvalues
    for j in range(1, k_subset + 1):
        s[j] = -1 - 2 * theta[j]
        
    eigenvalues_s = list(s.values())

    # Step 4: Find the LCM of the absolute values of the eigenvalues of S
    lcm_factors = Counter()
    for val in eigenvalues_s:
        val_factors = get_prime_factorization(val)
        for prime, exponent in val_factors.items():
            lcm_factors[prime] = max(lcm_factors[prime], exponent)

    print("The maximum order among all elements of the Smith group of S is the least common multiple of the absolute values of the eigenvalues of S.")
    
    # Calculate the terms of the final multiplication equation
    terms = []
    result = 1
    for prime, exponent in sorted(lcm_factors.items()):
        term = prime ** exponent
        terms.append(str(term))
        result *= term
        
    # Step 5: Print the final equation and the result
    equation_str = " * ".join(terms)
    print("\nThe equation for the maximum order is:")
    print(f"{equation_str} = {result}")

    return result

# Execute the solution
final_answer = solve_and_print()
# The final answer is wrapped according to the format requirements.
# <<<3083244199898525>>>