def solve_stress():
    """
    This function determines the stressed syllable in a list of Old Russian phrases
    based on a derived set of rules.
    """

    # List of phrases to analyze
    phrases = [
        'i ne znali',
        'i povelo že',
        'ne vymyla že',
        'ponesla',
        'vyvela že',
        'i unesli'
    ]

    # The parts of the phrases that are determined to be stressed based on the rules.
    stressed_parts = [
        'zna',  # Rule 5: Root stress
        'lo',   # Rule 3: End stress
        'vy',   # Rule 1: 'vy-' prefix stress
        'la',   # Rule 3 (extended): Feminine end stress
        'vy',   # Rule 1: 'vy-' prefix stress
        'li'    # Rule 3: End stress
    ]
    
    explanations = [
        "In 'i ne znali', the root 'zna' is stressed according to the root-stress pattern.",
        "In 'i povelo že', the verb follows the end-stress pattern, placing stress on the suffix 'lo'.",
        "In 'ne vymyla že', the prefix 'vy-' is always stressed, so it takes precedence.",
        "In 'ponesla', the feminine verb form '-la' receives the stress, following the end-stress pattern for this type.",
        "In 'vyvela že', the prefix 'vy-' is always stressed, taking precedence over other parts.",
        "In 'i unesli', the verb stem 'u-nes-' follows the end-stress pattern, placing stress on the suffix 'li'."
    ]

    final_digits = []

    print("Analysis of Old Russian Stress:")
    print("-" * 30)

    for i in range(len(phrases)):
        phrase = phrases[i]
        part = stressed_parts[i]
        explanation = explanations[i]

        vowels = "aeiouy"
        
        # Find the absolute starting position of the stressed part in the phrase.
        start_index_of_part = phrase.find(part)
        
        # Find the index of the first vowel within that stressed part.
        first_vowel_index_in_part = -1
        for j, char in enumerate(part):
            if char in vowels:
                first_vowel_index_in_part = j
                break
        
        # Calculate the absolute index of the stressed vowel in the entire phrase.
        abs_vowel_index = start_index_of_part + first_vowel_index_in_part
        
        # Count vowels from the beginning to find the syllable number.
        syllable_count = 0
        for k, char in enumerate(phrase):
            if char in vowels:
                syllable_count += 1
                if k == abs_vowel_index:
                    result = syllable_count
                    break
        
        final_digits.append(str(result))
        print(f"Phrase: '{phrase}'")
        print(f"Rule: {explanation}")
        print(f"The stressed syllable is number: {result}")
        print("-" * 30)

    final_answer = "".join(final_digits)
    print(f"The combined sequence of stressed syllable numbers is: {final_answer}")
    print(f"\n<<<342314>>>")

solve_stress()