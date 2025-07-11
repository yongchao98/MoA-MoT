import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    Returns 0 if k > n or k < 0.
    """
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

def get_lambda_i(i, n, k, t):
    """
    Calculates the eigenvalue lambda_i for the Johnson graph J(n,k,t).
    i: index of the eigenspace (from 0 to k)
    n: size of the base set
    k: size of the subsets
    t: required intersection size for adjacency
    """
    val = 0
    for j in range(min(i, t) + 1):
        term = combinations(i, j) * \
               combinations(k - i, t - j) * \
               combinations(n - k - i, k - t - i + j)
        if j % 2 == 1:
            val -= term
        else:
            val += term
    return val

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def main():
    """
    Main function to solve the problem.
    """
    n = 50
    k = 5
    t = 3

    # Calculate eigenvalues lambda_i of the adjacency matrix A
    lambdas = [get_lambda_i(i, n, k, t) for i in range(k + 1)]
    
    # Calculate eigenvalues theta_i of the Seidel matrix S
    thetas = []
    
    # Eigenvalue on V_0
    num_vertices = combinations(n, k)
    theta_0 = num_vertices - 1 - 2 * lambdas[0]
    thetas.append(theta_0)
    
    # Eigenvalues on V_i for i > 0
    for i in range(1, k + 1):
        theta_i = -1 - 2 * lambdas[i]
        thetas.append(theta_i)

    # We need the LCM of the absolute values of the eigenvalues
    abs_thetas = [abs(t) for t in thetas]

    print("The absolute values of the eigenvalues of the Seidel matrix S are:")
    for i, val in enumerate(abs_thetas):
        print(f"|theta_{i}| = {val}")

    # Calculate the LCM of all absolute theta values
    result_lcm = 1
    for val in abs_thetas:
        result_lcm = lcm(result_lcm, val)
        
    equation_values = ", ".join(map(str, abs_thetas))
    print(f"\nThe maximum order is the LCM of these numbers.")
    print(f"LCM({equation_values}) = {result_lcm}")
    
    print("\n<<<" + str(result_lcm) + ">>>")

if __name__ == "__main__":
    main()
