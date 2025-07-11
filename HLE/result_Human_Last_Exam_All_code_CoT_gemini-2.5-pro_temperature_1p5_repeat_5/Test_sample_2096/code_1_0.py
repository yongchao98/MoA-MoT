def solve_riddle():
    """
    This script solves the riddle by breaking down its clues.
    """

    # --- Step 1: Identify the clues from the riddle ---
    pope_name = "Paul II"
    key_decade = "1960s"
    
    # --- Step 2: Extract the numerical hint from the Pope's name ---
    # The Roman numeral "II" in the Pope's name is the first part of the puzzle.
    hint_number_1 = 2
    # The riddle uses this number to form the number 22.
    final_number = 22

    # --- Step 3: Connect the number to a work from the specified decade ---
    # We look for a famous work from the 1960s containing the number 22.
    # Joseph Heller's novel "Catch-22" was published in 1961.
    year_of_publication = 1961
    solution_name = "Catch-22"

    # --- Step 4: Print the logical steps and the final answer ---
    print("Solving the riddle step-by-step:")
    print(f"The first clue is the name of the Pope: {pope_name}.")
    print(f"The Roman numeral 'II' provides the number: {hint_number_1}")
    print(f"This number is used to form the key number in the solution: {final_number}")
    print(f"The second clue is that 'X' was written in the {key_decade}.")
    print(f"The famous work 'Catch-22' was published in {year_of_publication}, which fits the clue.")
    print("\nTherefore, 'X' is the term that describes a shameful, paradoxical situation.")
    print(f"The final answer is: {solution_name}")

solve_riddle()