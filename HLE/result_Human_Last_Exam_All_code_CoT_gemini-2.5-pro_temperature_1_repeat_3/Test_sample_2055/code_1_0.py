import math
from collections import defaultdict

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of num.
    """
    factors = defaultdict(int)
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] += 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] += 1
    return factors

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k".
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

def calculate_s_eigenvalues():
    """
    Calculates the eigenvalues of the Seidel matrix S.
    """
    v = 50
    k = 5
    d = 2 # Adjacency is for intersection size k-d=3

    # Calculate eigenvalues of the adjacency matrix A
    theta = []
    for j in range(k + 1):
        val = 0
        for i in range(d + 1):
            term = ((-1)**i * combinations(j, i) *
                    combinations(k - j, d - i) *
                    combinations(v - k - j, d - i))
            val += term
        theta.append(val)
        
    n = combinations(v, k)
    k_deg = theta[0]

    # Calculate eigenvalues of the Seidel matrix S
    sigma = []
    # Principal eigenvalue
    sigma.append(n - 1 - 2 * k_deg)
    # Other eigenvalues
    for j in range(1, k + 1):
        sigma.append(-1 - 2 * theta[j])
        
    return sigma

def solve():
    """
    Calculates the maximum order of an element in the Smith group of S.
    """
    # Step 1 & 2: Calculate eigenvalues of S
    sigma = calculate_s_eigenvalues()
    
    print("Eigenvalues of S:")
    for i, s_val in enumerate(sigma):
        print(f"sigma_{i} = {s_val}")
    print("-" * 20)

    # Step 3: Compute differences and find max p-adic valuations
    differences = []
    for i in range(len(sigma)):
        for j in range(i + 1, len(sigma)):
            diff = abs(sigma[i] - sigma[j])
            differences.append(diff)
            
    max_exponents = defaultdict(int)
    for diff in differences:
        if diff == 0:
            continue
        factors = get_prime_factorization(diff)
        for p, exp in factors.items():
            if exp > max_exponents[p]:
                max_exponents[p] = exp

    # Step 4: Calculate the final answer
    max_order = 1
    
    equation_parts = []
    for p, exp in sorted(max_exponents.items()):
        term = p ** exp
        max_order *= term
        equation_parts.append(f"{p}^{exp}")

    print("The maximum order is the product of the following prime powers:")
    print(" * ".join(equation_parts))
    print(f"\nResult: {max_order}")
    
    return max_order

solve()
