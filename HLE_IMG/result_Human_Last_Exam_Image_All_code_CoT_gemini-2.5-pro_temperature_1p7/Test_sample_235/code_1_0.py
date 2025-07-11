import re
try:
    import pronouncing
except ImportError:
    print("Please install the 'pronouncing' library to run this script: pip install pronouncing")
    exit()

def analyze_poem_meter():
    """
    Analyzes the metric pattern of the given erasure poem.
    """
    # 1. Define the poem's text from the image
    poem_text = "rules and lines, an intricate spider's web work"
    print(f"Analyzing the poem: \"{poem_text}\"")
    print("-" * 40)

    # 2. Process words to find syllables and stress patterns
    # The 'pronouncing' library requires word lookups, so we handle special cases.
    words = re.findall(r"[\w']+", poem_text.lower())
    
    # The CMU Pronouncing Dictionary doesn't have possessives like "spider's",
    # but "spiders" has the same syllabification and stress, so we map it.
    word_lookup_map = {
        "spider's": "spiders"
    }

    syllable_counts = []
    stress_pattern_list = []

    print("Syllable and Stress Analysis per Word:")
    for word in words:
        lookup_word = word_lookup_map.get(word, word)
        pronunciations = pronouncing.phones_for_word(lookup_word)
        
        if not pronunciations:
            print(f"  - '{word}': Could not find in dictionary. Skipping.")
            continue
        
        # Use the first pronunciation variant found
        first_pronunciation = pronunciations[0]
        count = pronouncing.syllable_count(first_pronunciation)
        stresses = pronouncing.stresses(first_pronunciation)
        
        print(f"  - '{word}': {count} syllable(s), stress pattern '{stresses}'")
        syllable_counts.append(count)
        stress_pattern_list.append(stresses)

    # 3. Calculate and display total syllables and the full stress pattern
    total_syllables = sum(syllable_counts)
    syllable_equation = " + ".join(map(str, syllable_counts))
    full_stress_pattern = "".join(stress_pattern_list)
    
    print("-" * 40)
    print("Overall Poetic Structure:")
    print(f"Total Syllable Calculation: {syllable_equation} = {total_syllables} syllables")
    print(f"Full Stress Pattern (1=stress, 0=no stress): {full_stress_pattern}")
    print("-" * 40)

    # 4. Evaluate against the provided metric pattern choices
    print("Evaluating Against Answer Choices:")
    
    # Choice A: Free Verse - The default if no other pattern fits.
    
    # Choice B: Iambic Pentameter (10 syllables, 5x 'da-DUM')
    is_iambic_pentameter = total_syllables == 10
    print(f"B. Iambic Pentameter (10 syllables): Poem has {total_syllables}. -> {'Match' if is_iambic_pentameter else 'No Match'}")

    # Choice C: Alexandrine (12 syllables)
    is_alexandrine = total_syllables == 12
    print(f"C. Alexandrine (12 syllables): Poem has {total_syllables}. -> {'Match' if is_alexandrine else 'No Match'}")

    # Choice D: Sapphic (11 syllables, specific meter)
    # The classic English Sapphic line has 11 syllables.
    # Meter: Trochee | Trochee | Dactyl | Trochee | Spondee
    # Stress:   10    |   10    |  100   |   10    |   11
    sapphic_pattern = "10101001011"
    is_sapphic = (total_syllables == 11 and full_stress_pattern == sapphic_pattern)
    print(f"D. Sapphic (11 syllables, complex meter): Poem has {total_syllables} syllables.")
    print(f"   The required pattern is {sapphic_pattern}. The poem's pattern is {full_stress_pattern}. -> {'Match' if is_sapphic else 'No Match'}")

    # Choice E: Loose Iambic Trimeter (~6 syllables)
    is_trimeter = 5 <= total_syllables <= 7
    print(f"E. Loose Iambic Trimeter (~6 syllables): Poem has {total_syllables}. -> {'Match' if is_trimeter else 'No Match'}")
    
    # Choice F: American Sentence (17 syllables)
    is_american_sentence = total_syllables == 17
    print(f"F. American Sentence (17 syllables): Poem has {total_syllables}. -> {'Match' if is_american_sentence else 'No Match'}")
    
    print("-" * 40)

    # 5. Conclude the analysis
    if is_sapphic:
        print("Conclusion: The poem's single line has 11 syllables and perfectly fits the stress pattern of a Sapphic line.")
    else:
        print("Conclusion: The poem does not fit a classical meter and would be considered free verse.")

if __name__ == '__main__':
    analyze_poem_meter()