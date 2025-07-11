def solve_bulgakov_puzzle():
    """
    This function solves a literary puzzle based on Bulgakov's 'The Master and Margarita'.

    1.  The Moscow character tormented in Chapter 18 from the provided list is Sokov,
        the barman who is told he will die of liver cancer. While a sparrow torments
        Poplavsky in the same chapter, Poplavsky is not on the list, making Sokov
        the only plausible choice from the options.
    2.  The corresponding bird in the Jerusalem chapters is a barn swallow.
    3.  This barn swallow flies around Pontius Pilate in Chapter 26, "The Burial",
        while he is suffering from a migraine and guilt.
    """
    
    moscow_character = "Sokov"
    jerusalem_bird = "barn swallow"
    jerusalem_character = "Pontius Pilate"

    # The final output combines the three elements as requested.
    # The print statement explicitly shows each part of the final answer.
    print(f"{moscow_character};{jerusalem_bird};{jerusalem_character}")

solve_bulgakov_puzzle()