import math

def C(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def gcd(a, b):
    """Computes the greatest common divisor of a and b"""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b"""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def calculate_max_order():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """
    n = 50
    k = 5

    # Step 1: Calculate eigenvalues of the Seidel matrix S.
    # The matrix A is the adjacency matrix of the graph where vertices are adjacent
    # if their intersection size is 3. This is the graph A_2 of the Johnson scheme J(50, 5).
    # S = J - I - 2*A_2
    
    # Calculate eigenvalues P_2j of A_2 on eigenspaces V_j
    P_2 = []
    for j in range(k + 1):
        p2j = C(k - j, 2) * C(n - k - j, 2) \
              - j * (k - j) * (n - k - j) \
              + C(j, 2)
        P_2.append(p2j)

    # The valency (degree) of the graph A_2 is P_2[0].
    valency = P_2[0]
    
    # Total number of vertices
    num_vertices = C(n, k)
    
    # Calculate eigenvalues theta_j of S = J - I - 2*A_2
    eigenvalues = []
    
    # For j=0, the eigenvector is the all-ones vector.
    # Eigenvalue of J is num_vertices, eigenvalue of I is 1, eigenvalue of A_2 is its valency.
    theta_0 = num_vertices - 1 - 2 * valency
    eigenvalues.append(theta_0)
    
    # For j > 0, the all-ones vector is orthogonal to the eigenspace.
    # So eigenvalue of J is 0.
    for j in range(1, k + 1):
        theta_j = -1 - 2 * P_2[j]
        eigenvalues.append(theta_j)
        
    print("The distinct eigenvalues of the Seidel matrix S are:")
    for i, val in enumerate(eigenvalues):
        print(f"theta_{i}: {val}")
    print("-" * 20)

    # Step 2: Compute absolute differences between all pairs of eigenvalues.
    differences = []
    print("The absolute differences |theta_i - theta_j| are:")
    for i in range(len(eigenvalues)):
        for j in range(i + 1, len(eigenvalues)):
            diff = abs(eigenvalues[i] - eigenvalues[j])
            differences.append(diff)
            print(f"|{eigenvalues[i]: >7} - ({eigenvalues[j]: >7})| = {diff}")
    print("-" * 20)
    
    # Step 3: Find the LCM of these differences.
    if not differences:
        max_order = 0
    else:
        result_lcm = differences[0]
        for i in range(1, len(differences)):
            result_lcm = lcm(result_lcm, differences[i])
        max_order = result_lcm
        
    print(f"The maximum order is the least common multiple (LCM) of these differences.")
    print(f"Final Answer: {max_order}")
    
calculate_max_order()
