import syllables

def analyze_poem_meter():
    """
    Analyzes the metric pattern of the erasure poem by counting syllables per line.
    """
    lines = [
        "rules",
        "and lines",
        "an intricate spider's web",
        "work"
    ]

    print("Analyzing the poem's structure by counting syllables per line:\n")
    line_syllable_counts = []
    for i, line_text in enumerate(lines):
        # Clean up the line for word splitting
        words = line_text.replace(',', '').split()
        line_total_syllables = 0
        
        # Calculate syllables for the whole line
        for word in words:
            line_total_syllables += syllables.estimate(word)
        
        line_syllable_counts.append(line_total_syllables)
        print(f"Line {i+1}: '{line_text}' -> {line_total_syllables} syllables")

    print(f"\nThe syllable pattern of the lines is: {line_syllable_counts}")
    
    print("\n--- Conclusion ---")
    print("The poem is structured into four lines with irregular lengths (1, 2, 7, and 1 syllables).")
    print("This lack of a consistent metrical pattern or uniform line length is the defining characteristic of free verse.")
    print("While the total syllable count of all words is 11 (the length of a Sapphic line), the poem's presentation with intentional, uneven line breaks points to free verse as the overarching form.")

analyze_poem_meter()