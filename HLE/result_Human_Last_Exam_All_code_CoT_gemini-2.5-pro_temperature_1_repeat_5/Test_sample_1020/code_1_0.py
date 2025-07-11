import math

def calculate_smallest_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: An integer.
        d: An odd integer such that 2 <= d <= n.
    """
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraints not met. Required: 2 <= d <= n. Got: n={n}, d={d}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be an odd integer. Got: d={d}.")
        return

    # The smallest complexity is given by the formula:
    # 2 + 2 * sum_{i=(d+1)/2}^{d-1} P(n, i)
    # where P(n, i) is the number of i-permutations of n.
    
    start_i = (d + 1) // 2
    end_i = d - 1
    
    perm_sum = 0
    perm_values = []
    
    for i in range(start_i, end_i + 1):
        # Calculate P(n, i) = n! / (n-i)!
        p_n_i = math.perm(n, i)
        perm_sum += p_n_i
        perm_values.append(p_n_i)

    total_complexity = 2 + 2 * perm_sum
    
    # Build the equation string for the final output
    sum_str = " + ".join(map(str, perm_values))
    if len(perm_values) > 1:
        sum_str = f"({sum_str})"
    elif len(perm_values) == 0: # This happens if d=1, but d>=2 is a constraint
        sum_str = "0"
    
    equation = f"2 + 2 * {sum_str} = {total_complexity}"
    
    print(f"For n={n} and d={d}, the smallest complexity is:")
    print(equation)

# Example usage with n=10 and d=5
# You can change these values to compute for other cases.
n_example = 10
d_example = 5
calculate_smallest_complexity(n_example, d_example)

# Another example with n=6 and d=3
n_example_2 = 6
d_example_2 = 3
calculate_smallest_complexity(n_example_2, d_example_2)
