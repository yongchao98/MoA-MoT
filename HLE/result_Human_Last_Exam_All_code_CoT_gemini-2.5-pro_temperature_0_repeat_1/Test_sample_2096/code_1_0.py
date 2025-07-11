def solve_riddle():
    """
    This function solves the riddle by breaking it down logically and printing the steps.
    """
    print("This is a riddle that relies on wordplay. Here is the step-by-step solution:")
    
    # Step 1: Identify the main clue from the Pope's name.
    pope_name = "Paul II"
    clue_from_name = "II"
    print(f"1. The central figure in the riddle is {pope_name}. The crucial part of his name is the Roman numeral '{clue_from_name}'.")

    # Step 2: Determine what would be a shameful condition for a Pope.
    shameful_condition = "the inability to read and write"
    print(f"2. The riddle states that being 'X' was considered shameful for a Pope. A very shameful condition for a learned figure like a Pope would be {shameful_condition}.")

    # Step 3: Identify the word for that condition.
    answer_word = "illiterate"
    print(f"3. The single word for this condition is '{answer_word}'.")

    # Step 4: Connect the word to the clue from the Pope's name.
    print(f"4. The word '{answer_word}' contains the letters 'i' twice, which is a pun on the Roman numeral '{clue_from_name}' from the Pope's name.")
    
    # Step 5: Point out the misdirection.
    print("5. The mention of the 1960s and the fall of Constantinople is a red herring designed to distract the solver.")

    print("\nTherefore, the final answer is:")
    print(answer_word)

solve_riddle()