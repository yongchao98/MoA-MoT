def solve_old_russian_stress():
    """
    This script determines the stressed syllable in Old Russian phrases based on a set of derived rules.
    It analyzes each phrase for specific components (prefixes, particles, roots) and applies a
    hierarchical rule system to find the correct syllable number.
    """

    phrases_to_analyze = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    # Verb roots and their inherent stress classes
    verb_classes = {
        'zna': 'A',  # Class A: Root-stressed
        'my': 'A',
        'nes': 'B',  # Class B: Ending-stressed
        've': 'B'
    }
    
    # Known prefixes (excluding the special case 'vy')
    prefixes = ['po', 'u']

    def get_syllable_number(phrase_nospace, component):
        """
        Calculates the syllable number of a component's vowel.
        A syllable is defined by a vowel (a, e, i, o, u, y). The number is the
        1-based count of vowels from the start of the phrase.
        """
        vowels = "aeiouy"
        
        # Find the starting character index of the component in the phrase
        if component.startswith('-'):  # For endings like '-la'
            search_str = component[1:]
            start_index = phrase_nospace.rfind(search_str)
        else:
            start_index = phrase_nospace.find(component)
        
        if start_index == -1: return "Error: component not found"

        # Find the absolute character index of the first vowel in or after the component start
        vowel_abs_pos = -1
        for i in range(start_index, len(phrase_nospace)):
            if phrase_nospace[i] in vowels:
                vowel_abs_pos = i
                break

        if vowel_abs_pos == -1: return "Error: vowel not found"

        # Count vowels up to and including the found vowel's position
        syllable_count = 0
        for i in range(vowel_abs_pos + 1):
            if phrase_nospace[i] in vowels:
                syllable_count += 1
        return syllable_count

    def get_stressed_syllable(phrase):
        """
        Applies the hierarchical rules to find the stressed syllable number.
        """
        phrase_nospace = phrase.replace(" ", "")

        # Rule 1: The prefix 'vy-' is always stressed.
        if 'vy' in phrase:
            return get_syllable_number(phrase_nospace, 'vy')

        # Rule 2: The particle 'ne' is present.
        if 'ne' in phrase:
            if phrase_nospace.endswith('la'): # Exception for feminine singular ending
                return get_syllable_number(phrase_nospace, '-la')
            else:
                return get_syllable_number(phrase_nospace, 'ne')

        # Rule 3: Default rules (no 'vy-' and no 'ne').
        root = None
        root_class = None
        for r, c in verb_classes.items():
            if r in phrase_nospace:
                root = r
                root_class = c
                break
        
        if root_class == 'A': # Class A verbs are root-stressed.
            return get_syllable_number(phrase_nospace, root)
        
        if root_class == 'B': # Class B stress depends on other elements.
            found_prefix = None
            for p in prefixes:
                # Ensure it's a prefix by checking its position
                if phrase_nospace.find(p) < phrase_nospace.find(root):
                    found_prefix = p
                    break
            
            if found_prefix: # If prefix is present, it's stressed.
                return get_syllable_number(phrase_nospace, found_prefix)
            else: # No prefix
                if 'že' in phrase: # If 'že' is present, it's stressed.
                    return get_syllable_number(phrase_nospace, 'že')
                else: # Otherwise, the ending is stressed.
                    # Find the last syllable.
                    vowels = "aeiouy"
                    last_vowel_pos = phrase_nospace.rfind(next(c for c in reversed(phrase_nospace) if c in vowels))
                    
                    syll_count = 0
                    for i in range(last_vowel_pos + 1):
                        if phrase_nospace[i] in vowels:
                             syll_count += 1
                    return syll_count
        
        return "Error: cannot determine stress"

    results = []
    print("Determining the stressed syllable for each phrase:")
    for phrase in phrases_to_analyze:
        stress_position = get_stressed_syllable(phrase)
        results.append(str(stress_position))
        print(f"'{phrase}' -> Stressed syllable is number: {stress_position}")
    
    final_answer = "".join(results)
    print("\nThe final sequence of digits is:")
    print(final_answer)

solve_old_russian_stress()