def calculate_feynman_diagrams(n):
    """
    Calculates a(n), which is the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n - 1)!!.

    This function prints the step-by-step calculation and returns the final value.
    """
    # The value for which to calculate the double factorial
    k = 2 * n - 1

    # Calculate the double factorial and create the equation string
    result = 1
    equation_terms = []
    
    # Iterate from k down to 1 with a step of -2
    for i in range(k, 0, -2):
        result *= i
        equation_terms.append(str(i))
    
    # Join the terms with ' * ' for the final output
    equation_str = " * ".join(equation_terms)
    
    # Print the full calculation process
    print(f"For n = {n}:")
    print(f"a({n}) = (2 * {n} - 1)!! = {k}!!")
    print(f"a({n}) = {k}!! = {equation_str} = {result}")

# The specific case requested is for a(3)
calculate_feynman_diagrams(3)