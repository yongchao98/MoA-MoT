def solve_riddle():
    """
    This script solves the riddle by breaking it down into logical steps.
    """
    # Step 1: Identify the historical figure and the 'shameful' quality.
    pope_1_name = "Pope Paul II"
    pope_1_century = "15th"
    accusation = "being an enemy of culture and learning"
    shameful_quality_adjective = "ILLITERATE"

    print("--- Step 1: Analyzing the historical subject ---")
    print(f"The subject is {pope_1_name}, who was Pope in the {pope_1_century} century.")
    print(f"His critics, the humanists, accused him of {accusation}.")
    print(f"This makes the shameful quality, 'X', the state of being: {shameful_quality_adjective}")
    print("\n")

    # Step 2: Analyze the 1960s clue via the name 'Paul'.
    pope_2_name = "Pope Paul VI"
    era = "1960s"
    document_year = 1967
    document_name = "Populorum progressio"
    key_word_in_document = "ILLITERACY"

    print("--- Step 2: Analyzing the 1960s clue ---")
    print(f"The clue mentions 'written in the {era}' and plays on the name 'Paul'.")
    print(f"This points not to Paul II, but to a different Pope: {pope_2_name}.")
    print(f"In the year {document_year}, {pope_2_name} wrote the encyclical '{document_name}'.")
    print(f"This document mentions the problem of: {key_word_in_document}")
    print("\n")

    # Step 3: Combine the clues to find the final one-word answer.
    print("--- Step 3: Final Conclusion ---")
    print(f"The 'shameful' accusation against Paul II was being '{shameful_quality_adjective}'.")
    print(f"The related word '{key_word_in_document}' was 'written' by Paul VI in {document_year}.")
    print("Therefore, the one-word answer 'X' is the adjectival form.")
    print("\n")
    print(f"Final Answer Equation:")
    print(f"Accusation against {pope_1_name} + Word written by {pope_2_name} in {era} = X")
    print(f"'{shameful_quality_adjective}' + '{key_word_in_document}' => {shameful_quality_adjective}")

solve_riddle()