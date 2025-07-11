def solve_riddle():
    """
    This script solves the riddle by identifying the crater, its name's etymology,
    and deriving the answer 'Y'.
    """

    # Step 1 & 2: The riddle points to Sagan Crater. The name to investigate is "Sagan".
    name = "Sagan"

    # Step 3: The etymology of the Hebrew name "Sagan" is "assistant" or "deputy",
    # specifically the title of the assistant to the High Priest.
    # This gives us the phrase "X of Y".
    etymological_meaning = "Assistant of the Priest"
    
    # Step 4: From the phrase, we can identify X and Y.
    X = "Assistant"
    Y = "the Priest"

    # Step 5: The puzzle asks for the value of Y.
    # The first part of the puzzle about what Carl Sagan wrote is a framing device.
    # The puzzle's solution lies in the literal etymology of his name.
    
    print(f"The name of the crater is '{name}'.")
    print(f"The etymological meaning of '{name}' fits the format 'X of Y'.")
    print(f"The phrase 'X of Y' is: '{etymological_meaning}'.")
    print(f"From this, we deduce:")
    print(f"X = '{X}'")
    print(f"Y = '{Y}'")
    print("\nThe answer to 'What is Y?' is:")
    print(Y)

solve_riddle()