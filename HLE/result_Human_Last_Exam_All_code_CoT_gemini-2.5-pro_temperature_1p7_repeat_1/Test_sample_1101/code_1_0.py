import math

def falling_factorial(x, n):
    """Computes the falling factorial (x)_n = x * (x-1) * ... * (x-n+1)"""
    if n == 0:
        return 1
    result = 1.0
    for i in range(n):
        result *= (x - i)
    return result

def binomial_coefficient(n, k):
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
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * (n choose m) * (n+1/2)_m
    and prints the expression.
    """
    if not isinstance(n, int) or n < 0:
        print("n must be a non-negative integer.")
        return

    total_sum = 0
    expression_parts = []

    for m in range(n + 1):
        # Calculate terms
        sign = (-1)**m
        binom_coeff = binomial_coefficient(n, m)
        
        # Calculate the falling factorial term: (n+1/2)_m
        # The ratio of Gamma functions is equivalent to this
        ff_val = falling_factorial(n + 0.5, m)
        
        term_val = sign * binom_coeff * ff_val
        total_sum += term_val

        # Prepare for printing
        sign_str = "-" if sign < 0 else "+"
        # For the first term, omit the plus sign
        if m == 0:
            sign_str = ""

        # Print the individual terms of the summation equation
        # Handle m=0 and m=1 cases for clarity
        if m == 0:
            expression_parts.append(f"{ff_val:.4f}")
        else:
            ff_parts = ' * '.join([f"({n+0.5 - i})" for i in range(m)])
            expression_parts.append(f" {sign_str} {binom_coeff} * {ff_val:.4f}")


    # Since the expression can be long, we don't print all multiplied components
    # Instead we show a simplified form.
    # The final equation would look like: 1.0000 - C(n,1)*term_1 + C(n,2)*term_2 - ... = total_sum
    final_expression = "".join(expression_parts)
    print(f"For n={n}:")
    
    full_equation_str = []
    sum_val = 0
    for m in range(n + 1):
        sign = (-1)**m
        binom = binomial_coefficient(n, m)
        # Using Gamma functions to show the original form
        term_gamma_val = math.gamma(n + 1.5) / math.gamma(n + 1.5 - m)
        sum_val += sign * binom * term_gamma_val
        
        op = "-" if m > 0 and sign > 0 else "+"
        if m==0 : op = ""
        elif sign < 0 and m > 0 : op = "-"

        #print component by component
        if m > 0:
           print(f" {op} ", end="")
        print(f"{binom}*{term_gamma_val:.2f}", end="")
    
    print(f" = {sum_val:.4f}")


# Example for a specific n, e.g., n=3
# You can change the value of n
n_to_solve = 4
solve_sum(n_to_solve)
