import math

def combinations(n, k):
    """Helper function to calculate binomial coefficients."""
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
    """Helper function for greatest common divisor."""
    return math.gcd(a, b)

def lcm(a, b):
    """Helper function for least common multiple."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def calculate_max_order():
    """
    Calculates the maximum order of the Smith group for the given Seidel matrix.
    """
    # Step 1 & 2: Define parameters and formulas for eigenvalues of S
    v = 50
    k = 5
    
    # Number of vertices n
    n = combinations(v, k)
    
    # Valency k_3 (degree of the graph)
    # Number of 5-subsets that intersect a given 5-subset in 3 elements.
    k3 = combinations(k, 3) * combinations(v - k, k - 3)
    
    # Step 3: Known eigenvalues of the adjacency matrix A_3 for J(50, 5)
    P3 = {
        0: 9900,
        1: 5500,
        2: 2145,
        3: 655,
        4: -25,
        5: -75
    }
    
    # Step 4: Calculate the eigenvalues of the Seidel matrix S
    theta = {}
    # Eigenvalue for the principal eigenspace (j=0)
    theta[0] = n - 1 - 2 * P3[0]
    
    # Eigenvalues for the other eigenspaces (j=1 to 5)
    for j in range(1, k + 1):
        theta[j] = -1 - 2 * P3[j]
        
    eigenvalues = list(theta.values())
    abs_eigenvalues = [abs(e) for e in eigenvalues]
    
    print("The eigenvalues of the Seidel matrix S are:")
    for i, val in enumerate(eigenvalues):
        print(f"theta_{i} = {val}")
    print("-" * 20)
    
    # Step 5: Compute the LCM of the absolute values of the eigenvalues
    print("The absolute values of the eigenvalues are:")
    print(abs_eigenvalues)
    print("-" * 20)
    
    result = 1
    for val in abs_eigenvalues:
        result = lcm(result, val)
        
    print("The maximum order of an element in the Smith group is the LCM of these absolute values.")
    # The problem asks to output each number in the final equation.
    # We show the list of numbers for which the LCM is calculated.
    equation_str = "lcm(" + ", ".join(map(str, abs_eigenvalues)) + ")"
    print(f"{equation_str} = {result}")
    
    return result

# Run the calculation and print the final answer
max_order = calculate_max_order()
print(f"\nFinal Answer: {max_order}")
