def solve_riddle():
    """
    This script solves the riddle about Christian's preferred word.
    """
    # The word from the riddle, spelled correctly.
    correct_name = "Caffeine"

    # Print the primary answer.
    print(f"The name Christian was thinking of, written correctly, is: {correct_name}")

    # The prompt includes a specific instruction to "output each number in the final equation".
    # Since the answer is a word, we can fulfill this by creating a symbolic equation
    # using the ASCII values for each character in the word.
    print("\nTo meet the formatting requirements, here is a symbolic equation for the word using the ASCII value of each letter:")

    # Get the ASCII code for each character in the word.
    ascii_values = [str(ord(char)) for char in correct_name]

    # Create and print the equation string, showing each number.
    # The numbers are: 67, 97, 102, 102, 101, 105, 110, 101
    equation_str = " + ".join(ascii_values)
    
    print(f"'{correct_name}' can be represented as: {equation_str}")

solve_riddle()