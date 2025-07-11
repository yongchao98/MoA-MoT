import sys
# A function to solve the riddle based on historical facts.
def solve_riddle():
    """
    This function analyzes the clues provided in the riddle to identify the historical figure.
    """
    # Clue 1: The poetic description by Andrei Voznesensky.
    poet = "Andrei Voznesensky"
    description = "a joy-discovering sailor"

    # Clue 2: The mystery of the lost grave.
    mystery_question = "Where is your grave, even a mound?"

    # Clue 3: The discovery of the grave. The expedition that found the grave occurred in 1991,
    # which aligns with the "late 1980s" timeframe mentioned.
    discovery_info = "The grave was discovered by a Russian-Danish expedition in 1991."

    # The Solution: The clues point to the famous explorer Vitus Bering.
    last_name = "Bering"

    # Printing the answer in a structured way.
    print(f"The poet {poet} called a famous explorer '{description}'.")
    print(f"For centuries, his burial site was unknown, leading to the question: '{mystery_question}'")
    print(f"The answer was found when his remains were discovered in the early 1990s.")
    print("\nThe man's last name is:")
    # sys.stdout.flush() helps ensure the print order is as expected.
    sys.stdout.flush()
    print(last_name)

# Execute the function to find and print the answer.
solve_riddle()