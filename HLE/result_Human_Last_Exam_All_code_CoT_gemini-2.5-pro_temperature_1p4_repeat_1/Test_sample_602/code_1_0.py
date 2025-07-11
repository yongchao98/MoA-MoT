import math

def get_l_of_n_expression():
    """
    Returns the symbolic expression for l(n) as a string.
    """
    expression = "2 + 2/n^2 - (2*(2*n - 1))/(n^2) * sqrt(n^2 - n + 1)"
    return expression

def main():
    """
    Prints the calculation steps and the final expression for l(n).
    """
    
    # The step-by-step derivation is provided in the text outside the code block.
    # Here, we will just print the final equation.

    final_expression_str = get_l_of_n_expression()
    
    # We print each number and symbol in the final equation.
    # Note: for variables and functions like 'n', 'sqrt', we print them as is.
    print("The exact value of l(n) is given by the equation:")
    print("l(n) = 2 + 2 / n^2 - (2 * (2 * n - 1)) / n^2 * sqrt(n^2 - n + 1)")

if __name__ == "__main__":
    main()
