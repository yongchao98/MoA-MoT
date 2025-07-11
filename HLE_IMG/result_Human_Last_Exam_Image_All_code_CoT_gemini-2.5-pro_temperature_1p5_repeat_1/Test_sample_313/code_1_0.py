def solve_riddle():
    """
    This script solves the riddle by calculating the numerical value
    of the identified phrase "Mad Hatter".
    """
    phrase = "MAD HATTER"
    # Remove space for calculation
    phrase_for_calc = phrase.replace(" ", "")

    numbers = []
    for char in phrase_for_calc:
        # Get the value for each letter (A=1, B=2, ...)
        # ord(char) - ord('A') gives a 0-based index, so we add 1.
        value = ord(char) - ord('A') + 1
        numbers.append(value)

    total = sum(numbers)

    # Create the equation string
    equation_str = " + ".join(map(str, numbers))

    print("The hidden phrase is 'Mad Hatter'.")
    print("Calculating its value by converting letters to numbers (A=1, B=2, etc.):")
    print(f"{phrase}:")
    print(f"{equation_str} = {total}")

solve_riddle()