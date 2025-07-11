def solve_latin_riddle():
    """
    This script solves the Latin riddle to determine the punishment for a thief
    caught in a garden by constructing a word from parts of four names.
    """

    # The riddle specifies the names from which to take the first syllables.
    # The names are Penelope, Dido, Cadmus, and Remus.
    part1_source = "Penelope"
    part1 = part1_source[:2]

    part2_source = "Dido"
    part2 = part2_source[:2]

    # For the pun to work, the riddle uses "Ca" from Cadmus.
    part3_source = "Cadmus"
    part3 = part3_source[:2]

    part4_source = "Remus"
    part4 = part4_source[:2]

    # The final word is the concatenation of these syllables.
    final_word = part1 + part2 + part3 + part4

    # The riddle asks for the word made "ex illis" (from those).
    # We will show the construction.
    print("The punishment is a word formed from four parts:")
    print(f"1. The first syllable of '{part1_source}': {part1}")
    print(f"2. The first syllable of '{part2_source}': {part2}")
    print(f"3. The first syllable of '{part3_source}': {part3}")
    print(f"4. The first syllable of '{part4_source}': {part4}")
    print("\nPutting the syllables together forms the equation for the punishment:")
    print(f"{part1} + {part2} + {part3} + {part4} = {final_word}")
    print(f"\nThe punishment is the Latin verb '{final_word}', which means 'to be sodomized'.")

solve_latin_riddle()