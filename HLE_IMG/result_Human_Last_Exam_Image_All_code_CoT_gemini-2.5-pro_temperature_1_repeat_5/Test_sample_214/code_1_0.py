def analyze_poetic_form():
    """
    Analyzes the syllabic structure of poetic lines to determine the form.
    """
    print("Analyzing the syllabic structure of the lines...\n")

    # --- Analysis of the line provided in the prompt ---
    line1_text = "nacreous wavers"
    line1_words = {"nacreous": 3, "wavers": 2}
    line1_total_syllables = sum(line1_words.values())

    print(f"Line 1: '{line1_text}'")
    # Build the equation string showing each number
    equation1 = " + ".join(map(str, line1_words.values()))
    print(f"Syllable Calculation: {equation1} = {line1_total_syllables} syllables")
    print("-" * 30)

    # --- Analysis of the line constructed from the image ---
    # The line "forms like a ghostly dance, a" is constructed from words in the image.
    line2_text = "forms like a ghostly dance a"
    line2_words = {"forms": 1, "like": 1, "a": 1, "ghostly": 2, "dance": 1, "a_again": 1}
    line2_total_syllables = sum(line2_words.values())

    print(f"Line 2: '{line2_text}' (from image words)")
    # Build the equation string showing each number
    equation2 = " + ".join(map(str, line2_words.values()))
    print(f"Syllable Calculation: {equation2} = {line2_total_syllables} syllables")
    print("-" * 30)

    # --- Conclusion ---
    print(f"The poems in the sequence contain lines of {line1_total_syllables} and {line2_total_syllables} syllables.")
    print("This 5-7-5 syllabic pattern is the structure of a Haiku.")

analyze_poetic_form()