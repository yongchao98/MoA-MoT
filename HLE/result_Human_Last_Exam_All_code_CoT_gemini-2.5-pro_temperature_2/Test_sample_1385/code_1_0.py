def solve_enclitic_order():
    """
    Determines and prints the fixed order of a given set of Old Russian enclitics.
    The order is based on established rules from historical Slavic linguistics.
    """

    # The given enclitics
    # бы (by) - conditional particle
    # же (zhe) - conjunctive/emphatic particle
    # еси (esi) - verbal clitic, "you are"
    # бо (bo) - conjunctive particle, "for, because"
    # мя (mya) - pronominal clitic, accusative "me"

    # The established hierarchical order:
    # 1. bo
    # 2. zhe
    # 3. by
    # 4. mya
    # 5. esi
    ordered_enclitics = ["бо", "же", "бы", "мя", "еси"]

    print("The correct order of the enclitics, forming a sequence, is:")

    # The user request asks to "output each number in the final equation".
    # Since we are dealing with words, we will output each word in the final sequence.
    print(ordered_enclitics[0], "-", ordered_enclitics[1], "-", ordered_enclitics[2], "-", ordered_enclitics[3], "-", ordered_enclitics[4])

solve_enclitic_order()