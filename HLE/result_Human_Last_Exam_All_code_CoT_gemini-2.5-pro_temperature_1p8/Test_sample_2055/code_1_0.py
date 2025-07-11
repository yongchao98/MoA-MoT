import math

def combinations(n, k):
    """Computes the number of combinations of k items from a set of n."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve():
    """
    Calculates the maximum order of the Smith group for the specified Seidel matrix.
    """
    # Step 1: Set parameters for the Johnson scheme J(v, k)
    v = 50
    k = 5

    # Step 2: Calculate eigenvalues of the adjacency matrix A of the graph Gamma
    A_eigenvalues = []
    for j in range(k + 1):
        # Formula for the eigenvalues of the A_2 matrix in the Johnson scheme J(v, k)
        term1 = combinations(k - j, 2) * combinations(v - k - j, 2)
        term2 = j * (k - j) * (v - k - j)
        term3 = combinations(j, 2)
        lambda_j = term1 - term2 + term3
        A_eigenvalues.append(lambda_j)

    # Step 3: Calculate eigenvalues of the Seidel matrix S
    num_vertices = combinations(v, k)
    lambda_0 = A_eigenvalues[0]

    S_eigenvalues = []
    # Eigenvalue for the principal eigenvector
    theta_0 = num_vertices - 1 - 2 * lambda_0
    S_eigenvalues.append(int(theta_0))

    # Eigenvalues for the other eigenspaces
    for j in range(1, k + 1):
        lambda_j = A_eigenvalues[j]
        theta_j = -1 - 2 * lambda_j
        S_eigenvalues.append(int(theta_j))

    # Step 4: The maximum order is the LCM of the absolute values of the eigenvalues.
    def gcd(a, b):
        return math.gcd(a, b)

    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b) if a != 0 and b != 0 else 0
    
    abs_S_eigenvalues = [abs(e) for e in S_eigenvalues]
    
    # Calculate the LCM of all these absolute values
    if not abs_S_eigenvalues:
        max_order = 0
    else:
        max_order = abs_S_eigenvalues[0]
        for i in range(1, len(abs_S_eigenvalues)):
            max_order = lcm(max_order, abs_S_eigenvalues[i])
            
    # Outputting the required information and the final answer.
    print("The distinct eigenvalues of the Seidel matrix S are:")
    print(", ".join(map(str, S_eigenvalues)))
    
    lcm_str = "lcm(" + ", ".join(map(str, abs_S_eigenvalues)) + ")"
    print(f"\nThe maximum order is computed as:\n{lcm_str}")
    
    print(f"\nResult: {max_order}")
    return max_order

final_answer = solve()
print(f"\n<<< {final_answer} >>>")