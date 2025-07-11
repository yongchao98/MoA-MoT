def solve_poetic_form():
    """
    Analyzes the syllable count of a poem's lines to determine its form.
    """
    # Line 4 analysis
    line4_text = "nacreous wavers"
    line4_words = {
        "nacreous": 3,
        "wavers": 2
    }
    line4_syllables = sum(line4_words.values())

    print(f"Analyzing Line 4: '{line4_text}'")
    # Build and print the equation for line 4
    line4_equation_parts = [f"{word}({syllables})" for word, syllables in line4_words.items()]
    line4_num_parts = [str(s) for s in line4_words.values()]
    print(f"Syllable count: {' + '.join(line4_equation_parts)} = {' + '.join(line4_num_parts)} = {line4_syllables} syllables.")
    print("-" * 20)

    # Line 3 analysis
    # From the image, we select the words that form a coherent 7-syllable line.
    line3_text = "ghostly velum like a dance"
    line3_words = {
        "ghostly": 2,
        "velum": 2,
        "like": 1,
        "a": 1,
        "dance": 1
    }
    line3_syllables = sum(line3_words.values())

    print(f"Analyzing Line 3: '{line3_text}'")
    # Build and print the equation for line 3
    line3_equation_parts = [f"{word}({syllables})" for word, syllables in line3_words.items()]
    line3_num_parts = [str(s) for s in line3_words.values()]
    print(f"Syllable count: {' + '.join(line3_equation_parts)} = {' + '.join(line3_num_parts)} = {line3_syllables} syllables.")
    print("-" * 20)

    # Conclusion
    print("Conclusion:")
    print(f"The poem's third line has {line3_syllables} syllables and its fourth line has {line4_syllables} syllables.")
    print("This 7-5 syllable pattern forms the last two lines of a classic 5-7-5 Haiku.")
    print("\nTherefore, the poetic form the text adheres to is a Haiku.")

solve_poetic_form()