def print_tiling_formula():
    """
    Prints the derived formula for the number of ways to tile the L-shape.
    The formula is F_n^2 + F_{n-1}^2.
    """
    # The components of the formula, as requested by the prompt to output each number.
    index1 = "n"
    power1 = 2
    index2 = "n-1"
    power2 = 2

    print("The number of ways to tile the shape for a given size 'n' is:")
    # Using ** for exponentiation, as is standard in Python.
    # The curly braces around index2 ensure 'n-1' is treated as a single block.
    print(f"F_{index1}**{power1} + F_{{{index2}}}**{power2}")
    print("\nWhere F_k represents the k-th Fibonacci number (with F_1=1, F_2=1).")

print_tiling_formula()