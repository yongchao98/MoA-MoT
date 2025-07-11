def f(n, m):
    """
    This function determines if the first player has a winning probability
    strictly more than 50% in the 2D-Generalized NIM game on a random
    n x m binary matrix.

    Args:
        n (int): The number of rows in the matrix.
        m (int): The number of columns in the matrix.

    Returns:
        int: 1 if the first player's winning probability is > 50%, 0 otherwise.
    """
    # Based on the analysis, the condition P(Win) > 0.5 is equivalent
    # to checking if the product of the dimensions is greater than 1.
    product = n * m
    result = product > 1
    
    # Convert boolean True/False to integer 1/0
    return 1 if result else 0

def analyze_and_print(n, m):
    """
    Analyzes a given n, m pair and prints the reasoning and result.
    """
    print(f"Analyzing for a {n} x {m} matrix...")
    
    # The final equation is the condition that determines the outcome.
    product = n * m
    result_bool = product > 1
    
    print(f"The condition for the first player to have a winning probability > 50% is n * m > 1.")
    print(f"For n={n} and m={m}, the equation is: {n} * {m} = {product}.")
    print(f"Is {product} > 1? {result_bool}.")
    
    final_value = f(n, m)
    print(f"Thus, the value of the function f({n}, {m}) is {final_value}.")
    print("-" * 30)

# Demonstrate with some examples
analyze_and_print(1, 1)
analyze_and_print(2, 3)
analyze_and_print(5, 1)
