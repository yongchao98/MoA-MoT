import math

def falling_factorial(x, n):
    """Computes the falling factorial (x)_n = x * (x-1) * ... * (x-n+1)"""
    if n == 0:
        return 1.0
    res = 1.0
    for i in range(n):
        res *= (x - i)
    return res

def combinations(n, k):
    """Computes the binomial coefficient n choose k"""
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

def solve_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n)/Gamma(1.5+n-m)
    and prints the equation.
    """
    total_sum = 0
    
    # The term Gamma(1.5+n)/Gamma(1.5+n-m) is equivalent to the falling factorial (n+0.5)_m
    x = n + 0.5
    
    equation_parts = []
    
    for m in range(n + 1):
        comb = combinations(n, m)
        term_val = falling_factorial(x, m)
        
        signed_term = ((-1)**m) * comb * term_val
        total_sum += signed_term
        
        # Format the output for each term in the equation
        sign = "-" if m % 2 != 0 else "+"
        # The first term's plus sign is handled separately
        if m == 0:
            sign = "" if signed_term < 0 else ""
            equation_parts.append(f"({comb})*({term_val:.4f})")
        else:
            sign_str = "-" if signed_term < 0 else "+"
            # A bit of logic to avoid double signs, e.g., "+ -5"
            if sign_str == "-":
                 equation_parts.append(f"- ({comb})*({term_val:.4f})")
            else:
                 equation_parts.append(f"+ ({comb})*({term_val:.4f})")


    print(f"For n={n}:")
    # The prompt asks to output each number in the final equation.
    # Let's print the sum structure.
    full_equation = " ".join(equation_parts)
    print(f"S_{n} = {full_equation} = {total_sum:.4f}")


# Example calculation for n=4
n_val = 4
solve_sum(n_val)
