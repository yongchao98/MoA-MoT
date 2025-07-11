def reveal_explorer_name():
    """
    This script solves the riddle by revealing the last name of the explorer.
    The clues point to Vitus Bering, a famous explorer in Russian service.
    His grave was lost for centuries and rediscovered by an expedition that began its
    work in the late 1980s, perfectly matching the riddle's description.
    This code will construct his name from the ASCII values of its letters.
    """
    # ASCII character codes for the last name B-E-R-I-N-G
    last_name_codes = [66, 69, 82, 73, 78, 71]

    print("Solving the equation to find the last name:")

    # We will print each number in the 'equation'
    for i, code in enumerate(last_name_codes):
        # Print the number
        print(code, end="")
        # Print a plus sign between numbers
        if i < len(last_name_codes) - 1:
            print(" + ", end="")

    # The result of the 'equation' is the name itself
    print(" = ", end="")

    # Now, we convert the codes back to characters to get the final answer
    result_name = ""
    for code in last_name_codes:
        result_name += chr(code)

    print(result_name)

reveal_explorer_name()