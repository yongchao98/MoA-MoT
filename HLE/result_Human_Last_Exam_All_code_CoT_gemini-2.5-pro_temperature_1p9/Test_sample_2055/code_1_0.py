import math

def p_i_j(i, j, n, k):
    """
    Calculates the eigenvalue of the matrix A_i on the j-th eigenspace
    of the Johnson association scheme J(n, k).
    """
    val = 0
    for u in range(i + 1):
        term = ((-1)**u * math.comb(j, u) * 
                math.comb(k - j, i - u) * 
                math.comb(n - k - j, i - u))
        val += term
    return val

def solve():
    """
    Solves the problem by calculating the eigenvalues of the Seidel matrix S
    and then finding the lcm of their absolute values.
    """
    n = 50
    k = 5
    # Adjacency is for intersection size 3, which is distance i = k - 3 = 2.
    adj_dist = 2
    
    # Calculate the total number of vertices
    N = math.comb(n, k)
    
    # Calculate eigenvalues p_2(j) of the adjacency matrix A
    A_eigenvalues = [p_i_j(adj_dist, j, n, k) for j in range(k + 1)]
    
    # Calculate eigenvalues theta_j of the Seidel matrix S
    S_eigenvalues = []
    # For j=0
    theta_0 = N - 1 - 2 * A_eigenvalues[0]
    S_eigenvalues.append(theta_0)
    # For j=1 to k
    for j in range(1, k + 1):
        theta_j = -1 - 2 * A_eigenvalues[j]
        S_eigenvalues.append(theta_j)

    # Get absolute values for lcm calculation
    abs_S_eigenvalues = [abs(val) for val in S_eigenvalues]
    
    # The maximum order in the Smith group is the lcm of the absolute eigenvalues
    # Python's math.lcm can handle multiple arguments since 3.9
    result = math.lcm(*abs_S_eigenvalues)
    
    # Print the equation as requested
    equation_str = f"lcm({', '.join(map(str, abs_S_eigenvalues))})"
    print(f"The calculation is: {equation_str}")
    print(f"The result is: {result}")

solve()