import math

def calculate_sum(n):
    """
    Calculates the sum S(n) = sum_{m=0 to n} (-1)^m * C(n, m) * Gamma(n+1.5)/Gamma(n-m+1.5)
    and prints the full equation.
    """

    def falling_factorial(x, k):
        """Calculates the falling factorial (x)_k."""
        if k == 0:
            return 1.0
        res = 1.0
        for i in range(k):
            res *= (x - i)
        return res

    total_sum = 0.0
    equation_parts = []

    print(f"To find the value of the sum for n={n}, we calculate each term:")
    print(f"S({n}) = sum_{{m=0}}^{n} (-1)^m * C({n},m) * (Gamma({n}+1.5) / Gamma({n}-m+1.5))")
    print("-" * 30)

    for m in range(n + 1):
        sign = (-1)**m
        comb = math.comb(n, m)
        
        # The ratio of Gamma functions is more accurately calculated as a falling factorial.
        # Gamma(n+1.5) / Gamma(n-m+1.5) = (n+0.5)_m
        ff_val = falling_factorial(n + 0.5, m)
        
        term_val = sign * comb * ff_val
        total_sum += term_val
        
        # Build the string for the equation
        sign_str = "-" if sign < 0 else "+"
        if m == 0:
            sign_str = ""
            
        part_str = f"{sign_str} {comb} * {ff_val:.4f}"
        equation_parts.append(part_str)

    # Print the full equation
    final_equation_str = " ".join(equation_parts).lstrip('+ ')
    print(f"S({n}) = {final_equation_str}")
    print(f"S({n}) = {total_sum:.4f}")


# Let's calculate for n = 4 as an example.
calculate_sum(4)
