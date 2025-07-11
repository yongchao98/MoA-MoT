import math

def combinations(n, k):
    """Helper function for combinations, returns 0 if k > n or k < 0."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def calculate_max_order():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """
    n_param = 50
    k_param = 5

    # Step 1: Calculate eigenvalues of the adjacency matrix A
    # The graph is the 2nd association graph of the Johnson scheme J(50, 5).
    # The eigenvalues lambda_j = p_2(j) are given by the formula:
    # p_2(j) = C(5-j,2)C(45-j,2) - j(5-j)(45-j) + C(j,2)
    # for j=0, ..., 5.
    lambda_eigenvalues = []
    for j in range(k_param + 1):
        term1 = combinations(k_param - j, 2) * combinations(n_param - k_param - j, 2)
        term2 = j * (k_param - j) * (n_param - k_param - j)
        term3 = combinations(j, 2)
        lambda_j = term1 - term2 + term3
        lambda_eigenvalues.append(lambda_j)

    # Step 2: Calculate eigenvalues of the Seidel matrix S = J - I - 2A
    num_vertices = combinations(n_param, k_param)
    
    seidel_eigenvalues = []
    # Principal eigenvalue
    eta_0 = num_vertices - 1 - 2 * lambda_eigenvalues[0]
    seidel_eigenvalues.append(eta_0)
    
    # Other eigenvalues
    for j in range(1, k_param + 1):
        eta_j = -1 - 2 * lambda_eigenvalues[j]
        seidel_eigenvalues.append(eta_j)
        
    abs_eigenvalues = [abs(e) for e in seidel_eigenvalues]

    # Step 3: Compute the LCM of the absolute values of the eigenvalues.
    # The maximum order in the Smith group is the largest invariant factor,
    # which is the LCM of the absolute values of the distinct eigenvalues.
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // math.gcd(a, b)

    result_lcm = 1
    for val in abs_eigenvalues:
        result_lcm = lcm(result_lcm, val)

    # Step 4: Print the results in the required format.
    print("The distinct eigenvalues of the Seidel matrix S are:")
    print(seidel_eigenvalues)
    print("\nThe maximum order of an element in the Smith group of S is the")
    print("least common multiple (LCM) of the absolute values of these eigenvalues.")
    
    equation_str = "lcm(" + ", ".join(map(str, abs_eigenvalues)) + ")"
    
    print(f"\n{equation_str} = {result_lcm}")

if __name__ == '__main__':
    calculate_max_order()
