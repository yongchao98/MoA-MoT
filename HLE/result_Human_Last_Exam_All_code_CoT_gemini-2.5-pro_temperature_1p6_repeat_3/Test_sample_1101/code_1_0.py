from fractions import Fraction
import math

def gamma_ratio(n, m):
    """
    Calculates Gamma(3/2 + n) / Gamma(3/2 + n - m) using Fractions for precision.
    This is the falling factorial [n + 1/2]_m.
    """
    if m == 0:
        return Fraction(1)
    
    # Start with n + 1/2
    val = Fraction(2 * n + 1, 2)
    result = Fraction(val)
    
    for i in range(1, m):
        val -= 1
        result *= val
    return result

def binomial_coeff(n, k):
    """
    Calculates the binomial coefficient C(n, k) using Fractions.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = Fraction(1)
    for i in range(k):
        res = res * (n - i) / (i + 1)
    return res

def calculate_s(n):
    """
    Calculates the sum S(n).
    """
    total_sum = Fraction(0)
    for m in range(n + 1):
        term_sign = (-1)**m
        binom = binomial_coeff(n, m)
        g_ratio = gamma_ratio(n, m)
        term = term_sign * binom * g_ratio
        total_sum += term
    return total_sum

def print_equation(n):
    """
    Prints the full equation for S(n).
    """
    print(f"For n={n}, the equation is:")
    
    equation_str = []
    
    for m in range(n + 1):
        sign = "-" if (-1)**m < 0 else "+"
        # first term does not need a sign unless it is negative
        if m == 0 and sign == '+':
            sign = ""
        
        binom_val = int(binomial_coeff(n, m))
        
        # Build the Gamma ratio part string
        gamma_expr = f"Gamma({n} + 3/2) / Gamma({n} + 3/2 - {m})"
        
        # Calculate the numeric value of the ratio for display
        ratio_val_frac = gamma_ratio(n, m)
        
        if ratio_val_frac.denominator == 1:
            ratio_val_str = str(ratio_val_frac.numerator)
        else:
            ratio_val_str = f"({ratio_val_frac})"
            
        term_str = f" {sign} C({n},{m}) * [{gamma_expr} = {ratio_val_str}]"
        equation_str.append(term_str)
        
    final_sum = calculate_s(n)
    
    print("S(n) =" + "".join(equation_str))
    print(f"\nS({n}) = {float(final_sum):.4f}")


# Let's compute for n=10 to demonstrate
n_to_calculate = 10
print_equation(n_to_calculate)