def get_stress_syllable(phrase):
    """
    Determines the stressed syllable number for a given Old Russian phrase
    based on a derived set of phonological rules.
    """
    words = phrase.split()
    verb_word = next((w for w in words if w not in ['i', 'ne', 'že']), None)

    def count_vowels(text):
        """Counts the number of vowels in a string to determine syllable count."""
        return sum(1 for char in text if char in "aeiouy")

    # Rule P: The prefix 'vy-' is always stressed.
    if verb_word.startswith('vy'):
        pre_vy_text = phrase[:phrase.find('vy')]
        return count_vowels(pre_vy_text) + 1

    # Rule N: Behavior with the 'ne' particle.
    if 'ne' in words:
        if verb_word.endswith('la'):  # Feminine ending
            # Stress on the final syllable of the verb.
            phrase_for_count = phrase.replace(' že', '')
            return count_vowels(phrase_for_count)
        else:  # Neuter or plural ending
            # Stress on the 'ne' particle.
            pre_ne_text = phrase[:phrase.find('ne')]
            return count_vowels(pre_ne_text) + 1

    # Rule D: Default stress patterns.
    is_root_stressed = any(root in verb_word for root in ['zna', 'my'])
    
    if is_root_stressed:
        # For root-stressed verbs, stress falls on the root.
        # Calculate syllables in any proclitics ('i') and prefixes ('po-') before the root.
        prefix = ''
        if verb_word.startswith('po'):
            prefix = 'po'
        
        pre_root_text = ''
        if 'i' in words:
            pre_root_text += 'i '
        pre_root_text += prefix
        
        return count_vowels(pre_root_text) + 1
    else:  # Ending-stressed verbs
        # Stress falls on the final syllable of the verb.
        phrase_for_count = phrase.replace(' že', '')
        return count_vowels(phrase_for_count)

def solve_and_print():
    """
    Analyzes a list of Old Russian phrases and prints the stressed syllable for each.
    """
    phrases_to_solve = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    results = []
    for p in phrases_to_solve:
        results.append(str(get_stress_syllable(p)))

    final_answer = "".join(results)
    print(final_answer)

solve_and_print()