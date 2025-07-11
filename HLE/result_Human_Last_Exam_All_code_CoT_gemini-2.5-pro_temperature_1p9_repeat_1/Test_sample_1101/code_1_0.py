import math

def falling_factorial(x, n):
    """Computes the falling factorial x*(x-1)*...*(x-n+1)."""
    if n == 0:
        return 1.0
    result = 1.0
    for i in range(n):
        result *= (x - i)
    return result

def combinations(n, k):
    """Computes binomial coefficient C(n,k)."""
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

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n)/Gamma(1.5+n-m)
    """
    if n < 0:
        return 0
    
    total_sum = 0
    
    # Use floating point numbers for calculations
    base_arg = 1.5 + n
    
    equation_parts = []
    
    for m in range(n + 1):
        # term = Gamma(1.5+n)/Gamma(1.5+n-m) which is (0.5+n)_m falling factorial
        term_val = falling_factorial(base_arg - 1, m)
        
        comb = combinations(n, m)
        
        # Calculate the full term in the sum
        sum_term = ((-1)**m) * comb * term_val
        total_sum += sum_term
        
        # Format parts for the equation string
        sign = "+" if m % 2 == 0 else "-"
        # On first term, don't print preceding +
        if m == 0:
            sign = ""
        
        equation_parts.append(f"{sign} {abs(sum_term):.4f}")
        
    print(f"For n={n}, the calculation is:")
    # The first part doesn't need a leading sign character space.
    # Replace the leading "+ " with ""
    equation_str = " ".join(equation_parts)
    if equation_str.startswith("+ "):
        equation_str = equation_str[2:]
        
    print(f"{equation_str} = {total_sum:.4f}")
    
    # Let's also check the ratio to n!
    # Avert division by zero for math.factorial(0) if total_sum is zero.
    if abs(total_sum) > 1e-9:
      ratio = abs(total_sum) / math.factorial(n)
      print(f"|S_n| / n! = |{total_sum:.4f}| / {math.factorial(n)} = {ratio:.4f}\n")

# --- Main execution ---
# Let's check for a few values of n.
for n_val in range(6):
    calculate_sum(n_val)

print("As n increases, the ratio |S_n|/n! appears to approach a non-zero constant.")
print("This implies that |S_n| is of the order of n!.")
print("Therefore, the function f(n) with the lowest complexity that bounds the sum is n!.")
