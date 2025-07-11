def solve_riddle():
    """
    This function solves the riddle by breaking it down into its logical components.
    """
    # Clue 1: A shameful trait/title for a Pope.
    shameful_role = "The head of a crime family"

    # Clue 2: A word 'written' in the 1960s. This points to a specific famous novel.
    year_of_writing = 1969
    famous_novel = "The Godfather"
    solution_word = "Godfather"

    print("Here is the step-by-step deduction to solve the riddle:")
    print("-" * 50)

    # Step 1: Identify the key clue connecting a shameful role to a specific decade.
    print(f"1. The riddle asks for a word 'X' that was 'written in the 1960s' and would be a 'shameful' title for a Pope.")

    # Step 2: Pinpoint the famous work from the 1960s.
    print(f"\n2. In the year {year_of_writing}, Mario Puzo wrote the famous novel '{famous_novel}'.")

    # Step 3: Explain why this word fits the 'shameful' description.
    print(f"\n3. The title of this novel, '{solution_word}', became associated with the leader of a criminal organization.")
    print(f"   This association would be considered shameful for a Pope, who is a spiritual leader.")

    # Step 4: Present the final logical equation.
    print("\n" + "=" * 50)
    print("The final logical equation is:")
    print(f"A word from a novel written in {year_of_writing} + A title shameful for a Pope = {solution_word}")
    print("=" * 50)

solve_riddle()
<<<Godfather>>>