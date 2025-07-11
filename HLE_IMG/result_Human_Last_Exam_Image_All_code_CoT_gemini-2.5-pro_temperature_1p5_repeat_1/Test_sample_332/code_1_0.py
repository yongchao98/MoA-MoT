import sys

def solve_riddle():
    """
    This script solves a riddle by deducing a phrase "X of Y" based on clues
    related to Carl Sagan, astronomy, and etymology.
    """

    # --- Clue Analysis ---
    # 1. Carl Sagan's concept: Images of Earth from space should inspire peace and unity,
    #    moving us away from war and conflict.
    # 2. A celestial feature's name has a meaning that fits the format "X of Y".

    # --- Hypothesis ---
    # The most fitting celestial name is "Mare Tranquillitatis" on the Moon.
    # Let's break down its name and meaning.
    latin_name = "Mare Tranquillitatis"
    
    # "Mare" is Latin for "Sea". This will be our X.
    X = "Sea"
    
    # "Tranquillitatis" is the genitive case of "tranquillitas", Latin for "tranquility".
    # It translates to "of Tranquility". This gives us Y.
    Y = "Tranquility"
    
    # --- Verification ---
    # The resulting phrase "Sea of Tranquility" perfectly encapsulates the sense of
    # peace and calm that Sagan hoped the cosmic perspective would bring to humanity,
    # serving as an ideal opposite to the "times of war" mentioned in the riddle.
    
    # The puzzle asks for the word 'Y'.
    final_answer = Y

    print(f"The riddle connects a Carl Sagan concept to the etymology of a crater's name.")
    print(f"The most logical fit is the lunar feature 'Mare Tranquillitatis'.")
    print(f"In Latin, 'Mare' means 'Sea' (X) and 'Tranquillitatis' means 'of Tranquility' (Y).")
    print(f"The full phrase, 'Sea of Tranquility', aligns with Sagan's message of peace.")
    print("-" * 20)
    print(f"The value of X is: {X}")
    print(f"The value of Y is: {Y}")
    print("-" * 20)
    print(f"The final answer is Y, which is '{final_answer}'.")

solve_riddle()