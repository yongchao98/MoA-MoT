import math

def combinations(n, k):
    """
    Calculates the binomial coefficient C(n, k).
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    # Using math.comb for efficiency and accuracy, requires Python 3.8+
    return math.comb(n, k)

def get_prime_factorization(num):
    """
    Returns the prime factorization of a number as a dictionary of {prime: exponent}.
    """
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def solve():
    """
    Main function to solve the problem.
    """
    # Step 1: Define graph parameters for the Johnson scheme J(v, k)
    v = 50
    k = 5
    # Adjacency is for intersection size 3, which corresponds to the matrix A_i with i = k - 3 = 2.
    i = 2

    print("The graph is from the Johnson scheme J(50, 5), with adjacency matrix A=A_2.")
    
    # Step 2: Calculate the eigenvalues of the adjacency matrix A
    # The eigenvalues P_i(j) of A_i are given by the formula:
    # P_i(j) = sum_{u=0 to i} (-1)^u * C(j,u) * C(k-j, i-u) * C(v-k-j, i-u)
    theta = []
    print("\nStep 1: Calculating eigenvalues of the adjacency matrix A...")
    for j in range(k + 1):
        val = 0
        for u in range(i + 1):
            term = ((-1)**u) * combinations(j, u) * combinations(k - j, i - u) * combinations(v - k - j, i - u)
            val += term
        theta.append(val)
        print(f"  theta_{j} = {val}")

    # Step 3: Calculate the eigenvalues of the Seidel matrix S = J - I - 2A
    n = combinations(v, k)
    seidel_eigenvalues = []
    
    print(f"\nStep 2: Calculating eigenvalues of the Seidel matrix S (n = {n})...")
    # Eigenvalue for the all-ones vector
    lambda_0 = n - 1 - 2 * theta[0]
    seidel_eigenvalues.append(lambda_0)
    print(f"  lambda_0 = n - 1 - 2*theta_0 = {n} - 1 - 2*{theta[0]} = {lambda_0}")

    # Eigenvalues for vectors orthogonal to the all-ones vector
    for j in range(1, k + 1):
        lambda_j = -1 - 2 * theta[j]
        seidel_eigenvalues.append(lambda_j)
        print(f"  lambda_{j} = -1 - 2*theta_{j} = -1 - 2*{theta[j]} = {lambda_j}")

    # Step 4: Calculate the maximum order, which is the lcm of the absolute values of the eigenvalues.
    # We do this by finding the highest power of each prime factor across all eigenvalues.
    print("\nStep 3: Calculating the maximum order in the Smith group of S.")
    print("This is the least common multiple (lcm) of the absolute values of the eigenvalues of S.")

    abs_eigenvalues = [abs(val) for val in seidel_eigenvalues if val != 0]

    max_prime_powers = {}
    for val in abs_eigenvalues:
        factors = get_prime_factorization(val)
        for p, a in factors.items():
            max_prime_powers[p] = max(max_prime_powers.get(p, 0), a)

    # Sort primes for a canonical representation of the final equation
    sorted_primes = sorted(max_prime_powers.keys())

    # Build the equation string and calculate the final result
    equation_parts = []
    max_order = 1
    for p in sorted_primes:
        a = max_prime_powers[p]
        if a > 1:
            equation_parts.append(f"{p}^{a}")
        else:
            equation_parts.append(f"{p}")
        max_order *= (p**a)
    
    equation_str = " * ".join(equation_parts)
    
    print("\nThe final calculation for the maximum order is:")
    print(f"  max_order = {equation_str}")
    print(f"  max_order = {max_order}")

solve()