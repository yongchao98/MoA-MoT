import math

def prime_factorization(n):
    """
    Returns a dictionary of prime factors of n.
    For example, prime_factorization(12) returns {2: 2, 3: 1}.
    """
    n = abs(n)
    if n < 2:
        return {}
    factors = {}
    # Check for factor 2
    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def main():
    """
    Main function to calculate the maximum order of the Smith group of S.
    """
    # Step 1: Define parameters for the Johnson scheme J(n, k)
    # Adjacency is based on intersection size k - i
    n = 50
    k = 5
    intersection_size = 3
    i = k - intersection_size  # This corresponds to the A_i matrix, so A_2

    # Step 2: Calculate the eigenvalues of the adjacency matrix A_i
    lambda_A = []
    print("Calculating eigenvalues of the adjacency matrix A_2...")
    for j in range(k + 1):
        val = 0
        for h in range(i + 1):
            # The general formula for eigenvalues of A_i in J(n,k)
            term = ((-1)**h) * math.comb(j, h) * math.comb(k - j, i - h) * math.comb(n - k - j, i - h)
            val += term
        lambda_A.append(val)
        print(f"  Eigenvalue for j={j}: lambda_{j}(A_2) = {val}")

    # Step 3: Calculate the eigenvalues of the Seidel matrix S = J - I - 2*A_2
    theta = []
    print("\nCalculating eigenvalues of the Seidel matrix S...")

    # Eigenvalue for the all-ones eigenvector (j=0)
    num_vertices = math.comb(n, k)
    valency = lambda_A[0]
    theta_0 = num_vertices - 1 - 2 * valency
    theta.append(theta_0)
    print(f"  Eigenvalue for j=0: theta_0 = {theta_0}")

    # Eigenvalues for other eigenspaces (j > 0)
    for j in range(1, k + 1):
        theta_j = -1 - 2 * lambda_A[j]
        theta.append(theta_j)
        print(f"  Eigenvalue for j={j}: theta_{j} = {theta_j}")
        
    # Step 4: Find the maximum p-adic valuation for each prime across all eigenvalues
    max_valuations = {}
    print("\nFinding prime factors of eigenvalues...")
    for t_val in theta:
        factors = prime_factorization(t_val)
        print(f"  Prime factorization of |{t_val}| is {factors}")
        for prime, power in factors.items():
            if prime not in max_valuations or power > max_valuations[prime]:
                max_valuations[prime] = power
    
    # Step 5: Calculate the maximum order, which is the product of p^max_valuation(p)
    print("\nThe maximum order is the product of the highest prime powers from the eigenvalues.")
    
    sorted_primes = sorted(max_valuations.keys())
    
    equation_parts = []
    for p in sorted_primes:
        equation_parts.append(f"{p}^{max_valuations[p]}")

    max_order = 1
    for p, power in max_valuations.items():
        max_order *= (p ** power)
    
    final_equation_str = " * ".join(equation_parts)
    print(f"The calculation is: {final_equation_str}")
    
    print(f"The maximum order is: {max_order}")


if __name__ == "__main__":
    main()
