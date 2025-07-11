def calculate_feynman_diagrams(n):
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n - 1)!!.
    """
    # Calculate the value k for the double factorial, which is 2n - 1
    k = 2 * n - 1

    # Initialize the result and a list to store the terms for the equation string
    result = 1
    equation_terms = []

    # Calculate the double factorial of k by iterating from k down to 1 with a step of -2
    for i in range(k, 0, -2):
        result *= i
        equation_terms.append(str(i))

    # Join the terms with " * " to create the equation string
    equation_str = " * ".join(equation_terms)

    # Print the explanation and the final calculated equation
    print(f"The number of diagrams a(n) for n={n} is given by (2*n - 1)!! = {k}!!")
    print(f"The calculation is:")
    print(f"{k}!! = {equation_str} = {result}")
    return result

# We need to find the value for a(3)
n_value = 3
a3 = calculate_feynman_diagrams(n_value)