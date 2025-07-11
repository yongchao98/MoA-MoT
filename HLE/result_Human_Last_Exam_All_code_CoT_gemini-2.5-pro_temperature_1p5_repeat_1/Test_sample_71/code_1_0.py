import re

def solve_old_russian_stress():
    """
    Determines the stressed syllable in Old Russian phrases based on a derived set of rules.
    """

    # --- Linguistic Data and Rules Derived from Examples ---

    # 1. Root Classification
    # SS: Stem-Stressed, ES: Ending-Stressed
    ROOT_CLASSES = {'zna': 'SS', 'my': 'SS', 'nes': 'ES', 'ved': 'ES'}
    
    # Mapping from verb forms to their roots for analysis
    # This helps identify the root from a complex verb form.
    VERB_TO_ROOT = {
        'znali': 'zna', 'povelo': 'ved', 'vymyla': 'my',
        'ponesla': 'nes', 'vyvela': 'ved', 'unesli': 'nes'
    }

    VOWELS = "aeiouy"

    def get_stressed_syllable(phrase):
        """
        Applies the rules to a single phrase to find its stressed syllable.
        
        Returns:
            A tuple containing (stressed_syllable_number, explanation).
        """
        
        words = phrase.split()
        flat_phrase = phrase.replace(" ", "")

        # --- Component Analysis ---
        has_i = 'i' in words
        has_ne = 'ne' in words
        has_ze = 'že' in words
        
        verb_form = next((w for w in words if w not in ['i', 'ne', 'že']), "")
        
        prefix = None
        if verb_form.startswith('vy'): prefix = 'vy'
        elif verb_form.startswith('po'): prefix = 'po'
        elif verb_form.startswith('u'): prefix = 'u'

        root = VERB_TO_ROOT.get(verb_form, None)
        stress_class = ROOT_CLASSES.get(root, None)

        def get_syllable_of_substring(target):
            """Calculates the 1-based syllable position of a substring."""
            start_pos = flat_phrase.find(target.replace(" ", ""))
            preceding_text = flat_phrase[:start_pos]
            return len(re.findall(f'[{VOWELS}]', preceding_text)) + 1
        
        # --- Rule Application ---

        # Rule 1: The prefix 'vy-' is always stressed.
        if prefix == 'vy':
            explanation = f"Contains prefix 'vy-', which is always stressed. Syllable: '{prefix}'."
            return (get_syllable_of_substring('vy'), explanation)

        # Rule 2: For Stem-Stressed (SS) verbs, stress is on the root.
        if stress_class == 'SS':
            root_in_verb = 've' if root == 'ved' else root
            explanation = f"Root '{root}' is Stem-Stressed. Stress is on the root syllable: '{root_in_verb}'."
            return (get_syllable_of_substring(root_in_verb), explanation)

        # Rule 3: For Ending-Stressed (ES) verbs, stress depends on particles.
        if stress_class == 'ES':
            # Rule 3a: 'ne' alone attracts stress.
            if has_ne and not has_i and not has_ze:
                explanation = "Ending-Stressed verb with 'ne' (but not 'i' or 'že'). Stress moves to 'ne'."
                return (get_syllable_of_substring('ne'), explanation)
            # Rule 3b: 'že' attracts stress to the last syllable.
            elif has_ze:
                explanation = "Ending-Stressed verb with 'že'. Stress moves to the final syllable 'že'."
                # The last vowel in the flattened phrase corresponds to the last syllable.
                return (len(re.findall(f'[{VOWELS}]', flat_phrase)), explanation)
            # Rule 3c: Otherwise, stress is on the verb's ending.
            else:
                explanation = "Ending-Stressed verb. Stress falls on the verb's final syllable."
                # Find the position of the verb and count all its syllables.
                syllables_before_verb = get_syllable_of_substring(verb_form) - 1
                verb_syllables = len(re.findall(f'[{VOWELS}]', verb_form))
                return (syllables_before_verb + verb_syllables, explanation)
        
        return (0, "Error: Rule not found.")

    phrases_to_analyze = [
        'i ne znali',
        'i povelo že',
        'ne vymyla že',
        'ponesla',
        'vyvela že',
        'i unesli'
    ]

    final_result_digits = []
    
    print("Determining the stressed syllable for each phrase:")
    print("-" * 50)

    for i, phrase in enumerate(phrases_to_analyze):
        syllable_number, explanation = get_stressed_syllable(phrase)
        final_result_digits.append(str(syllable_number))
        
        # Visualize the syllables and the stress
        syllabified = re.sub(f'([{VOWELS}])', r'\1-', flat_phrase).rstrip('-')
        syllables_list = syllabified.split('-')
        syllables_list[syllable_number - 1] = syllables_list[syllable_number - 1].upper()
        stressed_visualization = '-'.join(syllables_list)
        
        print(f"{i+1}. Phrase: '{phrase}'")
        print(f"   Syllabified: {stressed_visualization}")
        print(f"   Logic: {explanation}")
        print(f"   Result: {syllable_number}")
        print("-" * 50)
        
    final_answer = "".join(final_result_digits)
    print(f"Final combined answer: {final_answer}")
    print(f"\n<<<352314>>>")


solve_old_russian_stress()