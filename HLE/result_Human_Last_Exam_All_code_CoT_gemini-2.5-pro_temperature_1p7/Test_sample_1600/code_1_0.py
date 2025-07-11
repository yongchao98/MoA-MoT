def calculate_a(n):
    """
    Calculates a(n) = (2n-1)!!, the number of non-vanishing Feynman diagrams
    of order 2n.

    This function calculates the double factorial of (2n-1) and also
    prints the equation representing the calculation.
    """
    # The value to take the double factorial of
    k = 2 * n - 1

    # Initialize the result and the list to build the equation string
    result = 1
    equation_parts = []

    # Loop from k down to 1, with a step of -2 to get all odd numbers
    for i in range(k, 0, -2):
        result *= i
        equation_parts.append(str(i))

    # Join the numbers with " * " to form the multiplication string
    calculation_str = " * ".join(equation_parts)

    # Print the final equation with the result
    print(f"a({n}) = (2*{n}-1)!! = {k}!! = {calculation_str} = {result}")
    return result

# The problem asks for a(3)
n_value = 3
final_answer = calculate_a(n_value)