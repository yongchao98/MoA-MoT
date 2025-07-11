import math
from functools import reduce

def solve_seidel_smith_order():
    """
    Calculates the maximum order among all elements of the Smith group of the
    Seidel matrix S for the given graph.

    The graph Gamma has vertices as 5-element subsets of a 50-element set,
    with two subsets adjacent if their intersection has size 3.
    The Seidel matrix is S = J - I - 2A, where A is the adjacency matrix of Gamma.

    The maximum order of the Smith group is the largest invariant factor, d_N.
    For a matrix in the Bose-Mesner algebra of an association scheme with
    integer eigenvalues (like S), d_N is the LCM of the absolute values
    of the eigenvalues of S.
    """
    n = 50
    k = 5

    # Step 1: Calculate eigenvalues of the adjacency matrix A.
    # The graph described is A_2 in the Johnson scheme J(50,5).
    # Its eigenvalues are given by the formula:
    # P_2(j) = sum_{u=0 to 2} (-1)^u * C(j,u) * C(k-j,2-u) * C(n-k-j,2-u)
    
    def p_k(j, n, k, i):
        val = 0
        for u in range(i + 1):
            term = ((-1)**u) * math.comb(j, u) * math.comb(k - j, i - u) * math.comb(n - k - j, i - u)
            val += term
        return val

    # In our case, i=2 (intersection k-2 = 3).
    lambda_vals = [p_k(j, n, k, 2) for j in range(k + 1)]

    # Step 2: Calculate eigenvalues of the Seidel matrix S.
    # theta_0 corresponds to the all-ones eigenvector.
    num_vertices = math.comb(n, k)
    theta_0 = num_vertices - 1 - 2 * lambda_vals[0]

    # theta_j for j > 0 correspond to eigenvectors orthogonal to the all-ones vector.
    thetas = [theta_0]
    for j in range(1, k + 1):
        theta_j = -1 - 2 * lambda_vals[j]
        thetas.append(theta_j)

    # Step 3: Compute the LCM of the absolute values of the eigenvalues.
    abs_thetas = [abs(t) for t in thetas]

    def lcm(a, b):
        """Computes the least common multiple of a and b."""
        if a == 0 or b == 0:
            return 0
        return (a * b) // math.gcd(a, b)

    max_order = reduce(lcm, abs_thetas)
    
    # Final Output
    print("The eigenvalues of the Seidel matrix S are:")
    print(", ".join(map(str, thetas)))
    print("\nThe absolute values of these eigenvalues are:")
    numbers_str = ", ".join(map(str, abs_thetas))
    print(numbers_str)
    print("\nThe maximum order is the least common multiple (LCM) of these values.")
    print(f"LCM({numbers_str}) = {max_order}")
    print("\nThus, the final answer is:")
    print(max_order)


solve_seidel_smith_order()