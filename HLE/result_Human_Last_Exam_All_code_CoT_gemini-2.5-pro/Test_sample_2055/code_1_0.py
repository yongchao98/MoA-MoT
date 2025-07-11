import math

def combinations(n, k):
    """Computes the number of combinations 'n choose k'."""
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
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b) if a != 0 and b != 0 else 0

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary."""
    factors = {}
    d = 2
    temp_num = int(num)
    while d * d <= temp_num:
        while (temp_num % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
       factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def solve_graph_smith_order():
    """
    Calculates the maximum order among all elements of the Smith group of S.
    """
    n = 50
    k = 5

    # Step 1: Calculate eigenvalues of the adjacency matrix A
    A_eigenvalues = []
    for j in range(k + 1):
        val = (combinations(k, 2) * combinations(n - k - j, 0) * combinations(k - j, 0) -
               combinations(k - 1, 1) * combinations(n - k - j, 1) * combinations(k - j, 1) +
               combinations(k - 2, 2) * combinations(n - k - j, 2) * combinations(k - j, 2))
        
        # A simpler but equivalent formula for this specific graph (distance-2 Johnson graph)
        # P_{j,2} = C(5,2) - 4*C(5-j,1)*C(46-j,1) + C(5-j,2)*C(47-j,2)
        val_simple = 10 - 4 * (k - j) * (n - k - j) + combinations(k - j, 2) * combinations(n - k - 2 + j, 2)
        # After some algebra, the formula from the plan is correct:
        val_formula = 10 - 4 * (5 - j) * (46 - j) + combinations(5 - j, 2) * combinations(47 - j, 2)
        A_eigenvalues.append(val_formula)

    # Step 2: Calculate eigenvalues of the Seidel matrix S
    v = combinations(n, k)
    S_eigenvalues = []
    # For the all-ones eigenvector (j=0)
    s0 = v - 1 - 2 * A_eigenvalues[0]
    S_eigenvalues.append(s0)
    # For other eigenvectors (j>0)
    for j in range(1, k + 1):
        sj = -1 - 2 * A_eigenvalues[j]
        S_eigenvalues.append(sj)

    # Step 3: Compute the LCM of the absolute values of the eigenvalues of S
    abs_S_eigenvalues = [abs(s) for s in S_eigenvalues]
    
    result_lcm = 1
    for val in abs_S_eigenvalues:
        result_lcm = lcm(result_lcm, val)
        
    # Step 4: Format and print the result
    print(f"The distinct eigenvalues of the Seidel matrix S are: {S_eigenvalues}")
    
    print("\nThe maximum order is the least common multiple of the absolute values of these eigenvalues.")
    
    factors = get_prime_factorization(result_lcm)
    
    equation_parts = []
    for p, e in sorted(factors.items()):
        if e > 1:
            equation_parts.append(f"{p}^{e}")
        else:
            equation_parts.append(str(p))
            
    equation_str = " * ".join(equation_parts)
    
    print(f"\nThe final calculation is:\n{equation_str} = {result_lcm}")

solve_graph_smith_order()
<<<3083188729145425>>>