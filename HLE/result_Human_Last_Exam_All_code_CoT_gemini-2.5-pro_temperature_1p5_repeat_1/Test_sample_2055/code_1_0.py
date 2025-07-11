import math
from collections import defaultdict

def get_prime_factorization(num):
    """Computes the prime factorization of a positive integer."""
    factors = defaultdict(int)
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

def combinations(n, k):
    """Helper function for combinations, returns 0 if k > n or k < 0."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def main():
    """
    Solves the problem by calculating the eigenvalues of the Seidel matrix S
    and then finding the least common multiple of their absolute values,
    which corresponds to the largest invariant factor of S.
    """
    n = 50
    k = 5
    d = 2  # Adjacency for intersection size k-d=3, so d=2.

    # Calculate eigenvalues of the adjacency matrix A=A_2
    # Formula for eigenvalues P_d,i of A_d in J(n,k) is:
    # P_d,i = sum_{r=0 to d} (-1)^r * C(i,r) * C(k-i, d-r) * C(n-k-i, d-r)
    
    lambdas = []
    print("Step 1: Calculate eigenvalues of the adjacency matrix A.")
    for i in range(k + 1):
        lambda_i = 0
        for r in range(d + 1):
            term = ((-1)**r *
                    combinations(i, r) *
                    combinations(k - i, d - r) *
                    combinations(n - k - i, d - r))
            lambda_i += term
        lambdas.append(lambda_i)
        print(f"  λ_{i} = {lambda_i}")

    # Calculate eigenvalues of the Seidel matrix S = J - I - 2A
    num_vertices = combinations(n, k)
    
    thetas = []
    print("\nStep 2: Calculate eigenvalues of the Seidel matrix S.")
    # Eigenvalue for the all-ones eigenvector
    theta_0 = num_vertices - 1 - 2 * lambdas[0]
    thetas.append(theta_0)
    print(f"  θ_0 = num_vertices - 1 - 2*λ_0 = {num_vertices} - 1 - 2*{lambdas[0]} = {theta_0}")

    # Eigenvalues for eigenvectors orthogonal to the all-ones vector
    for i in range(1, k + 1):
        theta_i = -1 - 2 * lambdas[i]
        thetas.append(theta_i)
        print(f"  θ_{i} = -1 - 2*λ_{i} = -1 - 2*{lambdas[i]} = {theta_i}")

    # The maximum order in the Smith group is the largest invariant factor, d_n.
    # We approximate/assume this is the lcm of the absolute values of the eigenvalues.
    print("\nStep 3: Calculate the maximum order by finding the lcm of the absolute values of the eigenvalues of S.")
    
    # Dictionary to store the maximum exponent for each prime factor
    max_prime_powers = defaultdict(int)
    
    print("\nPrime factorization of each |θ_i|:")
    for theta in thetas:
        num = abs(theta)
        factors = get_prime_factorization(num)
        factor_str = " * ".join([f"{p}^{e}" for p, e in factors.items()])
        print(f"  |{theta}| = {num} = {factor_str}")
        for p, e in factors.items():
            if e > max_prime_powers[p]:
                max_prime_powers[p] = e

    print("\nMaximum powers of primes across all eigenvalues:")
    max_order = 1
    equation_parts = []
    # Sort primes for deterministic output
    sorted_primes = sorted(max_prime_powers.keys())
    for p in sorted_primes:
        power = max_prime_powers[p]
        print(f"  Prime {p}: max exponent is {power}")
        term = p**power
        max_order *= term
        equation_parts.append(f"{p}^{power}")

    final_equation = " * ".join(equation_parts)
    print(f"\nThe maximum order is the product of these highest prime powers:")
    print(f"Result = {final_equation}")
    print(f"Result = {max_order}")

if __name__ == "__main__":
    main()