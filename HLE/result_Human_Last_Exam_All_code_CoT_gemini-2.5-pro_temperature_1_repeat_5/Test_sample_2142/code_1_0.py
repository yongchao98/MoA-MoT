def solve_riddle():
    """
    Solves the Genesis P-Orridge riddle by identifying the words
    and printing the substitution cipher and the final answer.
    """
    # The two words are deduced from the theme of the riddle.
    Y = "GENDER"
    Z = "NONBINARY"

    # The full set of unique letters from both words.
    unique_letters = sorted(list(set(Y + Z)))

    # The "equation" is the substitution cipher mapping letters to numbers.
    # We map the sorted letters to digits 0-9.
    equation_map = {letter: i for i, letter in enumerate(unique_letters)}

    print("The final equation (the substitution cipher) is:")
    for letter, number in equation_map.items():
        print(f"{letter} = {number}")

    print("\nThe deciphered words are:")
    print(f"Y = {Y}")
    print(f"Z = {Z}")

    print("\nThe artist's commitment is to:")
    print(f"{Y} {Z}")

solve_riddle()
