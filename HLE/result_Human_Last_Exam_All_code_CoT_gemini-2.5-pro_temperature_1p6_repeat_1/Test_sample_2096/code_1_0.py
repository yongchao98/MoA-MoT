def solve_historical_riddle():
    """
    This script solves the riddle by breaking it down into logical steps
    and connecting the clues provided.
    """

    # --- Part 1: The core of the riddle ---
    print("Step 1: Analyzing the historical context of Pope Paul II.")
    pope_number = 2
    print(f"The riddle centers on Pope Paul {pope_number}.")
    print("Historically, Pope Paul II was criticized by humanist scholars as being hostile to learning.")
    print("Therefore, a shameful accusation against him would relate to a lack of education.")

    # --- Part 2: Identifying the shameful word 'X' ---
    print("\nStep 2: Determining the shameful word 'X'.")
    print("The most powerful insult for a Pope in that context would be to suggest he could not read or write.")
    shameful_word = "ILLITERATE"
    print(f"This makes the shameful word 'X': {shameful_word}.")
    print("The archbishop's statement is a clear defense against this specific charge.")

    # --- Part 3: Solving the final clue about the 1960s ---
    print("\nStep 3: Connecting the word to the 1960s.")
    decade_start = 1960
    print(f"The clue is that 'X' was 'written' in the {decade_start}s.")
    print("This refers to the global effort to combat illiteracy that was formally 'written' onto the world's agenda during that decade.")
    print("A key event was the 1965 UNESCO World Congress on the Eradication of Illiteracy.")

    # --- Part 4: The 'Final Equation' of Logic ---
    print("\nFinal Conclusion:")
    print("The logical 'equation' is as follows:")
    print(f"(Accusation shameful for Pope Paul {pope_number}) + (Concept 'written' onto the world stage in the {decade_start}s) = The Answer")
    print(f"'{shameful_word}' + (UNESCO's focus on Illiteracy in the 1960s) = {shameful_word}")
    print("\nThe word is therefore:")
    for char in shameful_word:
        print(char, end=' ')
    print() # for a new line

solve_historical_riddle()