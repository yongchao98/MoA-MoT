def solve_riddle():
    """
    This script solves a historical riddle by breaking it down into its components.
    """

    # --- Riddle Components ---
    pope = "Paul II"
    accusation = "Enemy of learning and culture"
    decade = "1960s"
    shameful_trait_concept = "Being unlearned during the Renaissance"

    # --- The logic connecting the components ---
    # The term for the accusation against the Pope, fitting the shame and 1960s clue.
    solution_word = "ILLITERATE"

    print("Solving the riddle step-by-step:")
    print("-" * 30)

    # The "equation" is the logical combination of the riddle's clues.
    # We will print each part of this logical equation.

    print(f"Part 1 (The Accusation): Pope {pope} was accused by humanists of being an '{accusation}'.")
    print(f"Part 2 (The Shame): This was considered shameful because a Pope was expected to be a patron of culture, not '{shameful_trait_concept}'.")
    print(f"Part 3 (The '1960s' Clue): The word for this trait was 'written' extensively in the {decade} during global campaigns, like UNESCO's World Congress on the Eradication of Illiteracy (1965).")

    print("-" * 30)
    print("The final equation is the sum of these logical parts:")
    print(f"'{accusation}' + '{shameful_trait_concept}' + 'Global keyword of the {decade}' = The final answer.")
    print("\nTherefore, the word 'X' is:")
    print(solution_word)


solve_riddle()
<<<ILLITERATE>>>