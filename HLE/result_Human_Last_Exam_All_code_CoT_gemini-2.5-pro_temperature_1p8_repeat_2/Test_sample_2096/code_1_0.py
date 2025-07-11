def solve_riddle():
    """
    This function solves the riddle by breaking down its components.
    """
    # Step 1: Analyze the historical figure.
    pope = "Paul II"
    pope_birth_name = "Pietro Barbo"

    # Step 2: Analyze the "shameful" trait based on his name.
    # The name "Barbo" comes from the Latin word "barba".
    meaning_of_barba = "beard"
    shameful_trait_inference = "bearded or hairy"
    reasoning = "In the 15th century, Latin Rite clergy were clean-shaven, in contrast to the bearded clergy of the Eastern Orthodox Church. So, being 'hairy' was seen as improper."

    # Step 3: Analyze the clue about the 1960s.
    # We need to find something "written" in the 1960s related to this trait.
    decade_clue = "1960s"
    famous_work_from_the_60s = "The rock musical 'Hair'"
    year_written = 1967

    # Step 4: Synthesize the clues to find the one-word answer.
    # The word "X" must describe the trait from Step 2 and relate to the work from Step 3.
    final_answer = "Hairy"

    print("Here is the step-by-step solution to the riddle:")
    print("-" * 50)
    print(f"1. The pope mentioned is {pope}, but his birth name is the real clue: {pope_birth_name}.")
    print(f"2. The name 'Barbo' is derived from the Latin word for '{meaning_of_barba}'.")
    print(f"3. Therefore, the trait considered 'shameful' for the Pope was being '{shameful_trait_inference}'.")
    print(f"4. The clue about a work 'written in the {decade_clue}' points to the famous rock musical '{famous_work_from_the_60s}', which premiered in {year_written}.")
    print(f"5. The word connecting the shameful trait ('bearded') and the 1960s musical ('Hair') is '{final_answer}'.")
    print("-" * 50)
    print(f"The final answer X is: {final_answer}")

solve_riddle()