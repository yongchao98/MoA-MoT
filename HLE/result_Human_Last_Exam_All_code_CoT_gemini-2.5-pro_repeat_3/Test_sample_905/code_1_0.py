def solve_monastery_riddle():
    """
    This function identifies a monastery based on historical clues
    and prints the answer.
    """
    # Clues from the user's question
    king = "King Philip II"
    year = 1585
    gate_name = "Golden Gate"
    reason_for_name = "the bronze plates covering it were gilded"
    insignias = "Sicily and Castile"

    # The identified monastery
    monastery = "Royal Monastery of Santa Maria of Guadalupe"

    # Constructing the final answer string, including the year from the prompt.
    answer = (
        f"The {monastery} has a {gate_name} where the insignias of {insignias} are displayed.\n"
        f"It was named by {king} because during his visit in {year}, {reason_for_name}."
    )

    print(answer)

solve_monastery_riddle()