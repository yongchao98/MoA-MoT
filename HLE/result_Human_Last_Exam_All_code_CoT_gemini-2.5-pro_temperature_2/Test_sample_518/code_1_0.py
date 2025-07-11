def solve_riddle():
    """
    Solves the Latin riddle by constructing a word from parts of names.
    """
    # The four names mentioned in the riddle.
    name1 = "Penelopes"
    name2 = "Didonis"
    name3 = "Cadmi"
    name4 = "Remi"

    # Extracting the first syllable (first two letters) from each name.
    # The riddle implies taking the first syllable from each.
    part1 = name1[:2]
    part2 = name2[:2]
    part3 = name3[:2]
    part4 = name4[:2]

    # Combine the parts to form the final word.
    punishment = part1 + part2 + part3 + part4

    print("The riddle builds the punishment word from four parts:")
    print(f"1. The first syllable of {name1}: {part1.upper()}")
    print(f"2. The first syllable of {name2}: {part2.upper()}")
    print(f"3. The first syllable of {name3}: {part3.upper()}")
    print(f"4. The first syllable of {name4}: {part4.upper()}")
    print("\nAssembling the final word, which is the punishment:")
    print(f"{part1.upper()} + {part2.upper()} + {part3.upper()} + {part4.upper()} = {punishment.upper()}")
    print(f"\nThe resulting Latin word is '{punishment.lower()}', a verb which means 'to sodomize'.")

solve_riddle()