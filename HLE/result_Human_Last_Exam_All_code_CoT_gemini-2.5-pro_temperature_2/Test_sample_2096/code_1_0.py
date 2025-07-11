def solve_riddle():
    """
    This function breaks down the riddle and prints the solution step-by-step.
    """
    # The riddle provides key numbers that help solve the puzzle.
    # We will use these numbers in our explanation.
    pope_number = 2  # From the name "Paul II"
    decade = 1960    # From the phrase "written in the 1960s"

    print("Analyzing the riddle by breaking it down into an 'equation' of clues:")
    print(f"Clue from number '{pope_number}': The mention of a Pope with this number in his name sets up a historical context.")
    print(f"Clue from number '{decade}': The answer, 'X', comes from a work written in this decade.")
    print("Clue from logic: It would be shameful for the 'Holy Father' (the Pope) to have the title 'X'.")
    
    print("\nSolving the equation:")
    print("1. A shameful title for a 'Father' figure would be the head of a crime family.")
    print(f"2. A famous novel about such a figure was written by Mario Puzo in 1969, which is in the {decade}s.")
    print("3. The title of this character is the solution.")

    answer = "Godfather"
    
    print("\n-------------------------")
    print(f"The word X is: {answer}")
    print("-------------------------")

solve_riddle()
<<<Godfather>>>