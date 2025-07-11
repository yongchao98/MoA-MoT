def solve_poetic_form_riddle():
    """
    Solves the riddle to determine the poetic form by combining clues
    from the prompt and the image to find a 5-7-5 syllable structure.
    """
    print("The poetic form can be determined by solving a riddle using numbers from the prompt and the image.")
    print("The target form is a Haiku, which has a 5-7-5 syllable structure.")
    print("-" * 20)

    # Deriving the first line's syllable count (5)
    num_from_image_left_1 = 2
    num_from_image_left_2 = 3
    first_line_syllables = num_from_image_left_1 + num_from_image_left_2
    print(f"1. The first line's syllable count comes from the numbers on the left margin of the image ({num_from_image_left_1} and {num_from_image_left_2}).")
    print(f"   Equation: {num_from_image_left_1} + {num_from_image_left_2} = {first_line_syllables} syllables.")
    print()

    # Deriving the second line's syllable count (7)
    num_from_prompt_1 = 3  # from "third line"
    num_from_prompt_2 = 4  # from "four-line stanza"
    second_line_syllables = num_from_prompt_1 + num_from_prompt_2
    print(f"2. The second line's syllable count comes from the numbers in the prompt's text ('third' -> {num_from_prompt_1}, 'four' -> {num_from_prompt_2}).")
    print(f"   Equation: {num_from_prompt_1} + {num_from_prompt_2} = {second_line_syllables} syllables.")
    print()

    # Deriving the third line's syllable count (5)
    syllables_nacreous = 3
    syllables_wavers = 2
    third_line_syllables = syllables_nacreous + syllables_wavers
    print(f"3. The third line's syllable count comes from the final poem 'nacreous wavers'.")
    print(f"   Syllables: 'nacreous' ({syllables_nacreous}) + 'wavers' ({syllables_wavers}).")
    print(f"   Equation: {syllables_nacreous} + {syllables_wavers} = {third_line_syllables} syllables.")
    print("-" * 20)

    print(f"The resulting syllable structure is {first_line_syllables}-{second_line_syllables}-{third_line_syllables}.")
    print("This 5-7-5 structure corresponds to the poetic form of a Haiku.")

solve_poetic_form_riddle()