import math

def calculate_feynman_diagrams(n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n - 1)!! (double factorial).
    This function also prints the full equation for the calculation.

    Args:
        n (int): The input value for the function a(n).
    """
    # The term for the double factorial is 2n - 1
    k = 2 * n - 1

    # Check for non-positive k which is not standard for this problem
    if k < 1:
        print(f"a({n}) = 0")
        return 0
    
    # Collect the terms for the equation string
    terms = []
    for i in range(k, 0, -2):
        terms.append(str(i))
        
    # Calculate the result using math.prod for Python 3.8+ or a loop
    # For compatibility, a loop is safer.
    result = 1
    for term in terms:
        result *= int(term)

    # Format the equation string
    equation_str = " * ".join(terms)
    
    # Print the final equation and result
    print(f"a({n}) = {equation_str} = {result}")

# The user is asking for a(3)
n_value = 3
calculate_feynman_diagrams(n_value)