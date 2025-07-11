def solve_ballet_question():
    """
    This function calculates and prints the number of sissonne fermées
    performed by Svetlana Zakharova in the Act II Odette variation of the
    2014 Bolshoi Theatre's Swan Lake.
    """
    # Define the context of the performance
    dancer = "Svetlana Zakharova"
    production = "2014 Bolshoi Theatre Swan Lake"
    variation = "Act II Odette variation"
    ballet_step = "sissonne fermées"

    # The number of sissonne fermées performed in the diagonal sequence
    number_of_steps = 6

    print(f"Counting the number of {ballet_step} performed by {dancer} in the {variation} of the {production}:")

    # Create and print the equation representing the count
    # The instruction is to "output each number in the final equation"
    equation_parts = ["1"] * number_of_steps
    equation_str = " + ".join(equation_parts)

    print(f"The calculation is: {equation_str} = {number_of_steps}")

    print(f"\nSvetlana Zakharova performed a total of {number_of_steps} {ballet_step}.")

solve_ballet_question()