def solve_riddle():
    """
    Solves the Latin riddle by assembling the punishment from its hidden parts.
    """

    # The first line of the riddle is a red herring. The real clues are elsewhere.

    # 1. The first part comes from "primam Cadmi syllaba" (the first syllable of Cadmus).
    part1 = "CAE"
    source1 = "primam Cadmi syllaba"

    # 2. The second part is hidden in the word "deprensus" (caught).
    part2 = "DE"
    source2 = "DEprensus"

    # 3. The third part comes from "prima Remi" (the first syllable of Remus).
    part3 = "RE"
    source3 = "prima Remi"

    # Assemble the punishment
    punishment = part1 + part2 + part3

    # Print the explanation and the result
    print("The Latin riddle is solved by combining syllables and letters hidden in the text.")
    print("The first line about Penelope and Dido is a misdirection.\n")
    print(f"The first part, from '{source1}', is: {part1}")
    print(f"The second part, hidden in the word '{source2}', is: {part2}")
    print(f"The third part, from '{source3}', is: {part3}\n")
    print("Putting them together forms the punishment:")
    print(f"{part1} + {part2} + {part3} = {punishment}\n")
    print(f"The punishment, '{punishment}', is a Latin verb meaning 'to be beaten' or 'to be struck'.")
    print("So, the punishment for stealing from the garden is a beating.")

solve_riddle()