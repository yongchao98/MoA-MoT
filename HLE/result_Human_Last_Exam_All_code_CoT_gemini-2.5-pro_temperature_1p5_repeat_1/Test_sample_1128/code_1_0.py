import math

def count_cool_strings(n):
    """
    Calculates the number of 'cool strings' of maximal length for n symbols.
    
    Args:
        n (int): The number of distinct symbols. Must be a positive integer.
        
    Returns:
        int: The total number of cool strings of maximal length 3n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return 0

    # The number of cool strings is given by the formula: (2n-1)! / (n-1)!
    # This can be computed as the product: n * (n+1) * ... * (2n-1)
    
    result = 1
    # Build the equation string for display
    equation_parts = []
    
    # The product runs from n to 2n-1
    for i in range(n, 2 * n):
        result *= i
        equation_parts.append(str(i))
        
    equation_str = " * ".join(equation_parts)
    
    print(f"For n = {n}, the number of cool strings is calculated as:")
    print(f"Calculation: {equation_str} = {result}")
    
    return result

# You can change the value of n here to test with other numbers.
n = 5
count_cool_strings(n)
