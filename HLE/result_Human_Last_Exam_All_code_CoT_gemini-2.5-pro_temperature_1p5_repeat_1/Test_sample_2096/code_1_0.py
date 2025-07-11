def solve_riddle():
    """
    Solves the riddle by identifying the pun and presenting the logic.
    """
    # The answer to the riddle.
    answer = "UNBELIEVER"

    # Explain the logic behind the riddle.
    print("This is a wordplay riddle. The key is not in the 15th-century history but in the clue about the 1960s.")
    print("1. The most shameful thing for a Pope is to be an 'UNBELIEVER'.")
    print("2. The clue is that 'X' was written in the 1960s.")
    print("3. This points to the 1966 hit song 'I'm a Believer' (by The Monkees).")
    print("4. The shameful word 'X' is the opposite of the song's premise.")
    print("-" * 20)
    print("The final answer is the word: " + answer)
    print("-" * 20)

    # To satisfy the instruction "output each number in the final equation",
    # we can represent the word as an 'equation' of its letter positions (A=1, B=2, ...).
    print("To represent this solution as a 'final equation' with numbers, we convert each letter to its position in the alphabet:")
    
    char_values = [ord(char.upper()) - ord('A') + 1 for char in answer]
    
    # Building the equation string with each number.
    equation_string = " + ".join(map(str, char_values))

    print(f"The equation for the word '{answer}' is:")
    # The final print statement outputs each number as requested.
    print(f"{equation_string}")

solve_riddle()