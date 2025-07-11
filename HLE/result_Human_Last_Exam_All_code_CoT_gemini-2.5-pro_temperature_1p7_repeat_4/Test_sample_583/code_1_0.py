import math

def f(n):
    """
    Calculates the value of f(n) based on the formula:
    f(n) = n + sum_{k=1}^{n-1} floor(log2(k))
    """
    if n <= 1:
        return 1
    
    # Calculate the sum of floor(log2(k)) for k from 1 to n-1
    sum_val = 0
    for k in range(1, n):
        sum_val += math.floor(math.log2(k))
        
    return n + sum_val

def calculate_limit_approximation(n):
    """
    Calculates the value of the expression f(n) / (n * log2(n)) for a given n.
    """
    if n <= 1:
        return float('nan')
        
    # Numerator of the expression
    numerator = f(n)
    
    # Denominator of the expression
    denominator = n * math.log2(n)
    
    # The final equation is ratio = numerator / denominator
    ratio = numerator / denominator
    
    print(f"Approximating the limit for n = {n}")
    print(f"The numerator f(n) is: {numerator}")
    print(f"The denominator n*log2(n) is: {denominator}")
    print(f"The ratio f(n)/(n*log2(n)) is: {ratio}")
    print("\nAs n approaches infinity, this ratio approaches 1.")


# Use a large value for n to approximate the limit
n_large = 100000
calculate_limit_approximation(n_large)

# The exact limit is 1. We can print the final answer as requested.
final_answer = 1
print(f"\nThe exact value of the limit is: {final_answer}")