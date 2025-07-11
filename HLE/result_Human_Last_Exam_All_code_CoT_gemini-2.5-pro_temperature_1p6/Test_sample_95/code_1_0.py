def solve_riddle():
    """
    This function solves the riddle by analyzing its clues.
    """

    # Clue 1: The setting is 19th-century industrial cities in Northern Europe.
    # The "smoky cities with restless air" refer to the heavy smog caused by the
    # Industrial Revolution, which would have drastically reduced long-distance visibility.
    deduction_1 = "Heavy smog in industrial cities obscured distant views."

    # Clue 2: A contrast is made with Milan.
    # This implies that "THEM" were visible from Milan. Milan is famously known
    # for its stunning view of a massive mountain range on clear days.
    deduction_2 = "The Alps are famously visible from Milan."

    # Clue 3: The name "Kasimir Graf" is likely a distraction.
    # The core of the riddle is the geographical and historical contrast.

    # Conclusion: The object(s) in question must be a large geographical feature
    # visible from Milan but easily obscured by smog. "THEM" is plural, which fits
    # the name of this mountain range.
    final_answer = "Alps"

    print(f"The answer is the mountain range visible from Milan but hidden by 19th-century industrial smog.")
    print(f"Therefore, 'THEM' are the: {final_answer}")

solve_riddle()