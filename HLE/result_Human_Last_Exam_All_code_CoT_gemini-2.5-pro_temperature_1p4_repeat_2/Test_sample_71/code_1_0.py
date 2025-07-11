def solve_stress():
    """
    This function determines the stressed syllable in Old Russian phrases based on a derived set of rules.
    It prints the analysis for each phrase and the final result.
    """
    
    # Phrases to analyze
    phrases = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    # Analysis based on the rules derived from the examples.
    # Class A roots (root stress): 'zna', 'my'
    # Class B roots (mobile stress): 'nes', 've' (from ved/vel)
    
    # Syllables (parts) for each phrase
    syllabified_phrases = [
        ['i', 'ne', 'zna', 'li'],                # i ne znali
        ['i', 'po', 've', 'lo', 'že'],           # i povelo že
        ['ne', 'vy', 'my', 'la', 'že'],          # ne vymyla že
        ['po', 'nes', 'la'],                     # ponesla
        ['vy', 've', 'la', 'že'],                # vyvela že
        ['i', 'u', 'nes', 'li']                  # i unesli
    ]
    
    # Rules determine the stressed part for each phrase
    stressed_parts = [
        'zna',   # Rule: 'zna' is Class A -> root stress.
        'po',    # Rule: 've' is Class B, prefixed, with 'že' -> prefix stress.
        'my',    # Rule: 'my' is Class A -> root stress.
        'la',    # Rule: 'nes' is Class B, prefixed, no particles -> ending stress (default).
        'vy',    # Rule: 've' is Class B, prefixed, with 'že' -> prefix stress.
        'li'     # Rule: 'nes' is Class B, prefixed, starts with 'i' -> ending stress.
    ]

    final_answer_digits = []
    
    print("Determining the stressed syllable for each phrase:")
    print("-" * 50)

    for i in range(len(phrases)):
        phrase = phrases[i]
        syllables = syllabified_phrases[i]
        stressed_syllable = stressed_parts[i]
        
        # Syllable number is the 1-based index of the stressed part
        stress_position = syllables.index(stressed_syllable) + 1
        final_answer_digits.append(str(stress_position))
        
        # Print the detailed analysis for each case
        syllable_str = "-".join(syllables)
        print(f"Phrase: '{phrase}'")
        print(f"Syllables: {syllable_str}")
        print(f"Stressed Syllable: '{stressed_syllable}'")
        print(f"Position: {stress_position}")
        print("-" * 50)
        
    final_answer = "".join(final_answer_digits)
    print(f"The final combined number is: {final_answer}")

solve_stress()
<<<323314>>>