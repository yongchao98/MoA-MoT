import math

def calculate_a(n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n-1)!!.
    """
    if n <= 0:
        return "n must be a positive integer."

    # Calculate the base for the double factorial
    base = 2 * n - 1

    # Calculate the double factorial using a loop
    result = 1
    terms = []
    current_term = base
    while current_term > 0:
        result *= current_term
        terms.append(str(current_term))
        current_term -= 2

    # Format the equation string showing all the numbers
    equation_str = " * ".join(terms)

    # Print the final result in the desired format
    print(f"a({n}) = ({2*n-1})!! = {equation_str} = {result}")
    return result

# Calculate and print a(3)
final_answer = calculate_a(3)