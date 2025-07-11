def solve_poetic_form_puzzle():
    """
    This function determines the poetic form by analyzing syllable counts
    and prints the solution.
    """
    # A Tanka is a five-line poem with a specific syllable structure.
    line1_syllables = 5
    line2_syllables = 7
    line3_syllables = 5
    line4_syllables = 7
    line5_syllables = 7

    # Calculate the total number of syllables
    total_syllables = (line1_syllables + line2_syllables + line3_syllables +
                       line4_syllables + line5_syllables)

    # Print the explanation and the final equation as requested.
    print("The poetic form is determined by its syllable structure.")
    print("The form consists of five lines.")
    print("The syllable structure is 5-7-5-7-7.")
    print("\nThe final equation representing the total syllables is:")
    print(f"{line1_syllables} + {line2_syllables} + {line3_syllables} + {line4_syllables} + {line5_syllables} = {total_syllables}")
    print("\nThis poetic form is known as a Tanka.")

solve_poetic_form_puzzle()