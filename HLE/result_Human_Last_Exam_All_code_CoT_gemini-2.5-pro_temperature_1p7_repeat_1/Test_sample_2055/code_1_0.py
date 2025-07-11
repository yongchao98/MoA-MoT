import math

def combinations(n, k):
    """
    Helper function to calculate combinations C(n, k).
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def get_lambda_j(j, n, k, i):
    """
    Calculates the eigenvalue lambda_j of the association matrix A_i
    of the Johnson scheme J(n, k).
    """
    val = 0
    for u in range(i + 1):
        term = (combinations(j, u) * 
                combinations(k - j, i - u) * 
                combinations(n - k - j, i - u))
        if u % 2 == 1: # for (-1)^u
            val -= term
        else:
            val += term
    return val

def lcm(a, b):
    """
    Computes the least common multiple of two integers.
    """
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b) if a != 0 and b != 0 else 0

def solve():
    """
    Main function to solve the problem.
    """
    # Parameters of the Johnson Scheme
    n = 50
    k = 5
    # Adjacency is for intersection size 3, which means we use matrix A_i
    # where k - i = 3 => 5 - i = 3 => i = 2
    i = 2

    print("Step 1: Calculate the eigenvalues of the adjacency matrix A.")
    lambdas = [get_lambda_j(j, n, k, i) for j in range(k + 1)]
    print(f"The eigenvalues of A are: {lambdas}")
    print("-" * 30)

    print("Step 2: Calculate the eigenvalues of the Seidel matrix S.")
    num_vertices = combinations(n, k)
    
    thetas = []
    # For j=0, the eigenvector is the all-ones vector
    theta_0 = num_vertices - 1 - 2 * lambdas[0]
    thetas.append(theta_0)
    
    # For j > 0, the eigenvectors are orthogonal to the all-ones vector
    for j in range(1, k + 1):
        theta_j = -1 - 2 * lambdas[j]
        thetas.append(theta_j)
    
    print(f"The eigenvalues of S are: {thetas}")
    abs_thetas = [abs(t) for t in thetas]
    print(f"The absolute values are: {abs_thetas}")
    print("-" * 30)

    print("Step 3: Calculate the maximum order in the Smith group.")
    print("This is the LCM of the absolute values of the eigenvalues of S.")
    
    result = 1
    for val in abs_thetas:
        result = lcm(result, val)

    # Output the final equation as requested
    lcm_expr = "lcm(" + ", ".join(map(str, abs_thetas)) + ")"
    print("\nFinal calculation:")
    print(f"{lcm_expr} = {result}")

    print("\nThus, the maximum order among all elements of the Smith group of S is:")
    print(result)
    
solve()