import math

def combinations(n, k):
    """Helper function for combinations, returns 0 if k > n or k < 0."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def compute_lambda_j(j, n, k, d):
    """
    Computes the eigenvalue lambda_j for the adjacency matrix A_d of the
    Johnson graph J(n,k). Adjacency in A_d means intersection size is k-d.
    """
    val = 0
    for i in range(d + 1):
        term = ((-1)**i * combinations(j, i) *
                combinations(k - j, d - i) *
                combinations(n - k - j, d - i))
        val += term
    return val

def solve():
    """
    Solves the problem by calculating eigenvalues and their LCM.
    """
    n = 50
    k = 5
    # Adjacency is for intersection size 3, which is k-2. So, d=2.
    d = 2

    # Calculate eigenvalues of the adjacency matrix A=A_2
    lambdas = [compute_lambda_j(j, n, k, d) for j in range(k + 1)]

    # Calculate eigenvalues of the Seidel matrix S = J - I - 2A
    thetas = []
    
    # Eigenvalue for the all-ones vector (j=0)
    N = combinations(n, k)
    k_d = lambdas[0] # valency of the graph
    theta_0 = N - 1 - 2 * k_d
    thetas.append(theta_0)

    # Eigenvalues for vectors orthogonal to the all-ones vector (j > 0)
    for j in range(1, k + 1):
        theta_j = -1 - 2 * lambdas[j]
        thetas.append(theta_j)

    print(f"The parameters of the Johnson graph are n={n}, k={k}.")
    print(f"The adjacency in the graph Gamma is defined by intersection size {k-d}.")
    print("\nThe eigenvalues of the adjacency matrix A of Gamma are:")
    for j, l in enumerate(lambdas):
        print(f"lambda_{j}: {l}")

    print("\nThe eigenvalues of the Seidel matrix S are:")
    for j, t in enumerate(thetas):
        print(f"theta_{j}: {t}")

    # Calculate the LCM of the absolute values of the eigenvalues of S
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // math.gcd(a, b)

    result_lcm = 1
    for t in thetas:
        # We need to compute lcm of integers, so convert to int
        result_lcm = lcm(result_lcm, int(t))

    print(f"\nThe maximum order among all elements of the Smith group of S is the")
    print("least common multiple of the absolute values of the eigenvalues of S.")
    print(f"\nEquation: max_order = lcm(|{thetas[0]}|, |{thetas[1]}|, |{thetas[2]}|, |{thetas[3]}|, |{thetas[4]}|, |{thetas[5]}|)")
    print(f"Result: {result_lcm}")
    
    # Return the final answer in the specified format
    print(f"\n<<<{result_lcm}>>>")

solve()