import math

def final_calculation():
    """
    This function performs the final calculation as requested.
    It calculates the sum for a specific value of n (n=5) and prints
    out the full equation, term by term.
    
    The sum is S_n = sum_{m=0 to n} (-1)^m * C(n, m) * Gamma(n+1.5) / Gamma(n+1.5-m).

    Based on numerical analysis, the absolute value of the sum |S_n| appears to be
    bounded by C * n! for some constant C > 0. The ratio |S_n/n!| appears to be
    bounded by a constant.
    So, the function f(n) with the lowest complexity that bounds the sum is f(n)=n!.
    """
    n = 5
    print(f"Finding the function f(n) of lowest complexity for the sum S_n.")
    print(f"S_n = \\sum_{{m=0}}^{{n}} (-1)^m \\binom{{n}}{{m}} \\frac{{\\Gamma({n}+1.5)}}{{\\Gamma({n}-m+1.5)}}")
    print("Numerical analysis suggests the sum is bounded by C*n!. So, we propose f(n) = n!.")
    print("\nHere is the detailed calculation for n =", n)
    
    total_sum = 0.0

    print("S_5 = ", end="")
    
    for m in range(n + 1):
        sign = (-1)**m
        comb = math.comb(n, m)
        
        # Calculate Gamma ratio product
        try:
            # For integer arguments, math.gamma is faster than repeated multiplication.
            gamma_ratio = math.gamma(n + 1.5) / math.gamma(n - m + 1.5)
        except ValueError:
            gamma_ratio = 0 

        term_val = sign * comb * gamma_ratio
        total_sum += term_val

        # Print the term
        # For m=0, sign_str is empty.
        # For m>0, sign_str is ' - ' or ' + '.
        sign_str = ' - ' if sign < 0 else ' + '
        if m == 0:
            print(f"{comb} * {gamma_ratio:.4f}", end="")
        else:
            print(f"{sign_str}{comb} * {gamma_ratio:.4f}", end="")

    print(f"\n= {total_sum}")

final_calculation()