def find_monastery():
    """
    This function stores and prints the answer to a historical trivia question
    based on a set of specific clues.
    """
    # Clues provided in the user's query
    clue_gate = "Has a Golden Gate."
    clue_insignias = "Displays the insignias of Sicily and Castile."
    clue_custom = "Kings would dismount before entering."
    clue_naming_monarch = "King Philip II"
    clue_naming_year = 1585
    clue_naming_reason = "The bronze plates on the gate were gilded for his visit."

    # The monastery that matches all the clues
    answer = "Royal Monastery of Santa María de Guadalupe"

    print("Searching for a monastery based on the following clues:")
    print(f"- {clue_gate}")
    print(f"- {clue_insignias}")
    print(f"- {clue_custom}")
    print(f"- Named 'Golden Gate' by {clue_naming_monarch} in {clue_naming_year}.")
    print(f"- Reason: {clue_naming_reason}")
    print("\n" + "="*40 + "\n")
    print("Answer:")
    print(f"The monastery is the {answer}.")

find_monastery()