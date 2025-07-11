import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    a = abs(a)
    b = abs(b)
    if a == 0 or b == 0:
        return 0
    return (a * b) // gcd(a, b)

def solve():
    """
    Calculates the maximum order of an element of the Smith group of the Seidel matrix S.
    """
    # Step 1: Define graph parameters
    n_set = 50
    k_subset = 5
    intersection_size = 3

    # Step 2: Calculate the eigenvalues of the adjacency matrix A
    A_eigenvalues = []
    print("The distinct eigenvalues of the adjacency matrix A are:")
    for j in range(k_subset + 1):
        theta_j = 0
        for i in range(j + 1):
            term = ((-1)**i * combinations(j, i) *
                    combinations(k_subset - i, intersection_size) *
                    combinations(n_set - k_subset - j + i, k_subset - intersection_size))
            theta_j += term
        A_eigenvalues.append(theta_j)
        print(f"theta_{j} = {theta_j}")

    # Step 3: Calculate the eigenvalues of the Seidel matrix S
    num_vertices = combinations(n_set, k_subset)
    
    S_eigenvalues = []
    # For j=0, the eigenvector is the all-ones vector
    sigma_0 = num_vertices - 1 - 2 * A_eigenvalues[0]
    S_eigenvalues.append(sigma_0)
    
    # For j > 0, the eigenvectors are orthogonal to the all-ones vector
    for j in range(1, k_subset + 1):
        sigma_j = -1 - 2 * A_eigenvalues[j]
        S_eigenvalues.append(sigma_j)

    print("\nThe distinct eigenvalues of the Seidel matrix S are:")
    for i, val in enumerate(S_eigenvalues):
        print(f"sigma_{i} = {val}")

    # Step 4: Find the maximum order in the Smith group by calculating the lcm
    # of the absolute values of the eigenvalues of S.
    max_order = 1
    for s_eig in S_eigenvalues:
        max_order = lcm(max_order, s_eig)
        
    print("\nThe final equation is the calculation of the least common multiple of the absolute values of these eigenvalues.")
    numbers_str = ", ".join(map(str, [abs(s) for s in S_eigenvalues]))
    print(f"max_order = lcm({numbers_str})")
    
    print(f"\nThe maximum order among all elements of the Smith group of S is: {max_order}")

solve()