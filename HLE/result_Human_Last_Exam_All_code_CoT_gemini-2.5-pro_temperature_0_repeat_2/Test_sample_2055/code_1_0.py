import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def main():
    """
    Calculates the maximum order of the Smith group for the specified Seidel matrix.
    """
    # Step 1: Define parameters
    n = 50
    k = 5
    d = 2  # Adjacency for intersection size k-d = 5-2 = 3

    # Step 2: Calculate eigenvalues of the adjacency matrix A = A_2
    eigenvalues_A = []
    for j in range(k + 1):
        # Formula for eigenvalues of A_d in J(n,k)
        # E_d(j) = sum_{i=0 to d} (-1)^i * C(j,i) * C(k-j, d-i) * C(n-k-j, d-i)
        term1 = combinations(j, 0) * combinations(k - j, d - 0) * combinations(n - k - j, d - 0)
        term2 = -combinations(j, 1) * combinations(k - j, d - 1) * combinations(n - k - j, d - 1)
        term3 = combinations(j, 2) * combinations(k - j, d - 2) * combinations(n - k - j, d - 2)
        E_d_j = term1 + term2 + term3
        eigenvalues_A.append(E_d_j)

    # Step 3: Calculate eigenvalues of the Seidel matrix S
    N = combinations(n, k)
    eigenvalues_S = []
    
    # Principal eigenvalue
    mu_0 = N - 1 - 2 * eigenvalues_A[0]
    eigenvalues_S.append(int(mu_0))
    
    # Other eigenvalues
    for j in range(1, k + 1):
        mu_j = -1 - 2 * eigenvalues_A[j]
        eigenvalues_S.append(int(mu_j))

    # Step 4 & 5: Compute the lcm of the absolute values of the eigenvalues of S
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // math.gcd(a, b)

    abs_eigenvalues_S = [abs(e) for e in eigenvalues_S]

    result_lcm = 1
    for val in abs_eigenvalues_S:
        result_lcm = lcm(result_lcm, val)

    print(f"The eigenvalues of the adjacency matrix A are: {eigenvalues_A}")
    print(f"The eigenvalues of the Seidel matrix S are: {eigenvalues_S}")
    
    # Format the final equation string
    equation_str = f"lcm({', '.join(map(str, abs_eigenvalues_S))})"
    
    print(f"The maximum order is the least common multiple of the absolute values of the eigenvalues of S.")
    print(f"{equation_str} = {result_lcm}")

if __name__ == "__main__":
    main()