def solve_riddle():
    """
    This function logically deciphers the provided riddle.
    """
    # 1. Analyze the geographical clues
    location_north = "smoky cities in the north of Europe"
    location_south = "Milan"
    print(f"Clue 1: Something is invisible from '{location_north}'.")
    print(f"Clue 2: The same thing IS visible from '{location_south}'.")
    print("Logic: This points to a very large geographical feature separating Northern Europe from Italy, where Milan is located.")
    print("-" * 30)

    # 2. Identify the feature based on geography
    feature = "The Alps"
    print(f"Analysis: The most famous mountain range visible from Milan but not from, for example, Germany, is {feature}.")
    print("On clear days, the Alpine mountain range is a well-known sight from the city of Milan.")
    print("-" * 30)

    # 3. Decipher the misdirection in the final clue
    person = "the German Kasimir Graf"
    reason = "he simply did not have enough imagination for this"
    print(f"Clue 3: {person} claims his inability to see THEM is because '{reason}'.")
    print("Logic: This is the riddle's clever misdirection. You do not need 'imagination' to see a mountain range; you need clear skies and proximity.")
    print("A German, living north of the Alps and in a smoggy industrial city, could not physically see them. He ironically blames his own imagination for this physical impossibility.")
    print("-" * 30)

    # 4. State the final conclusion
    print("Conclusion: The clues point to the great mountain range that stands as a barrier between the north and south of Europe.")
    print("Therefore, 'THEM' are...")
    print(feature)


solve_riddle()