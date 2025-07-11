import math

def gcd_list(numbers):
    """Computes the greatest common divisor of a list of numbers."""
    if not numbers:
        return 0
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = math.gcd(result, numbers[i])
    return result

def lcm(a, b):
    """Computes the least common multiple of two numbers."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b)

def lcm_list(numbers):
    """Computes the least common multiple of a list of numbers."""
    if not numbers:
        return 1
    result = numbers[0]
    for i in range(1, len(numbers)):
        result = lcm(result, numbers[i])
    return result

def solve():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """
    # Parameters of the Johnson scheme J(v, k)
    v = 50
    k = 5

    # Step 1: Calculate eigenvalues of the Johnson graph J(50,5) (distance-1 graph)
    theta = []
    print("Eigenvalues of the Johnson graph J(50,5):")
    for j in range(k + 1):
        val = (k - j) * (v - k - j) - j
        theta.append(val)
        print(f"  theta_{j} = {val}")

    # Step 2: Calculate eigenvalues of the adjacency matrix A of our graph (distance-2 graph)
    # The relation is A_2 = (1/4) * (A_1^2 - 48*A_1 - 225*I)
    p2 = []
    print("\nEigenvalues of the graph Gamma:")
    for j, t in enumerate(theta):
        val = (t**2 - 48 * t - 225) // 4
        p2.append(val)
        print(f"  p2_{j} = {val}")

    # Step 3: Calculate eigenvalues of the Seidel matrix S
    n = math.comb(v, k)
    deg_A = p2[0] # The degree of the graph is the first eigenvalue
    s0 = n - 1 - 2 * deg_A
    s_rest = [-1 - 2 * p for p in p2[1:]]
    s_eigenvalues = [s0] + s_rest
    print("\nEigenvalues of the Seidel matrix S:")
    for i, s in enumerate(s_eigenvalues):
        print(f"  s_{i} = {s}")

    # Step 4: Calculate the maximum order
    s_eigenvalues.sort(reverse=True)
    print("\nOrdered eigenvalues of S:")
    print(f"  {s_eigenvalues}")
    
    # Differences of consecutive ordered eigenvalues
    diffs = [s_eigenvalues[i] - s_eigenvalues[i + 1] for i in range(len(s_eigenvalues) - 1)]
    print("\nDifferences of consecutive ordered eigenvalues:")
    print(f"  {diffs}")

    # GCD of these differences
    g = gcd_list(diffs)
    print(f"\nGCD of differences (g) = {g}")
    
    # Terms for the final LCM calculation
    terms_for_lcm = [d // g for d in diffs]
    print("\nTerms for LCM calculation (differences / g):")
    print(f"  {terms_for_lcm}")

    # The maximum order is the LCM of these terms
    max_order = lcm_list(terms_for_lcm)
    
    # Final equation and result
    print(f"\nThe final equation is lcm{tuple(terms_for_lcm)}.")
    print(f"The maximum order is {max_order}.")

solve()