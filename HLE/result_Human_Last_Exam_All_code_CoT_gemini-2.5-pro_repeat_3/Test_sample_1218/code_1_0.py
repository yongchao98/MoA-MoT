import sympy

def calculate_max_n(k):
    """
    Calculates the maximum value of n in terms of k based on the derived formula.
    The formula is n = k^2 - k + 1.
    """
    if not isinstance(k, int) or k < 2:
        print("Error: k must be an integer greater than or equal to 2.")
        return None, None
    
    n = k**2 - k + 1
    equation = f"{k}*{k} - {k} + 1"
    return n, equation

def main():
    """
    Main function to demonstrate the calculation and print the result.
    """
    # The question is for a general k, so we can represent it symbolically.
    k_sym = sympy.Symbol('k')
    n_expr = k_sym**2 - k_sym + 1
    
    print("The maximum value of n in terms of k is given by the expression:")
    print(f"n = {n_expr}")
    print("\nLet's demonstrate this for some specific integer values of k:")
    
    test_cases = [2, 3, 4, 10]
    for k_val in test_cases:
        result, equation = calculate_max_n(k_val)
        if result is not None:
            # This formatting addresses the prompt's instruction to "output each number in the final equation"
            print(f"For k = {k_val}, the maximum value of n is {equation} = {result}")

if __name__ == "__main__":
    main()
