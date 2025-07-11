def calculate_feynman_diagrams(n):
    """
    Calculates the number of non-vanishing Feynman diagrams of order 2n.
    The formula is a(n) = (2n - 1)!!.
    This function prints the calculation step-by-step for a(n).
    """
    # Calculate the base of the double factorial
    k = 2 * n - 1

    print(f"The number of diagrams, a({n}), is given by the formula (2*n - 1)!!.")
    print(f"For n={n}, this is ({2}*{n} - 1)!! = {k}!!.")
    print("The calculation is:")

    # Calculate the double factorial and create the equation string
    result = 1
    equation_parts = []
    
    # Iterate from k down to 1 with a step of -2
    for i in range(k, 0, -2):
        result *= i
        equation_parts.append(str(i))

    equation_str = " * ".join(equation_parts)
    print(f"{equation_str} = {result}")
    
    return result

# We need to find a(3)
n_value = 3
final_answer = calculate_feynman_diagrams(n_value)