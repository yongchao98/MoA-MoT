def solve_history_question():
    """
    This function identifies the key figures and dates related to the
    stylization change of the French monarchy and outputs the answer.
    """
    # The monarch in question is Philip II.
    king_regnal_number = 2

    # The shift to "King of France" (Rex Franciae) began around 1190.
    year_of_title_change_start = 1190

    # The monarch's reign ended with his death in 1223, a year paired with the correct
    # biographer, Rigord, in the answer choices.
    year_of_kings_death = 1223

    # The chronicler who wrote the biography and gave Philip II his epithet "Augustus" was Rigord.
    biographer = "Rigord"

    # The problem asks for the year and the biographer. Choice C is (1223, Rigord).
    # We can create a simple equation to fulfill the prompt's requirements,
    # calculating the duration of the king's reign after the title change began.
    duration_after_change = year_of_kings_death - year_of_title_change_start

    print(f"The monarch was Philip II.")
    print(f"The biographer who provided his epithet was {biographer}.")
    print("The stylistic change to emphasize territoriality began around 1190.")
    print(f"The king's reign ended in {year_of_kings_death}.")
    print("\nBased on the options, the correct choice pairs the correct biographer with a key year of the reign.")
    print("We can represent the time from the title change to the end of the reign with the equation:")
    print(f"{year_of_kings_death} - {year_of_title_change_start} = {duration_after_change}")

solve_history_question()