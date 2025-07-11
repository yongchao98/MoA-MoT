def solve_old_russian_stress():
    """
    Determines the stressed syllable in Old Russian phrases based on a set of derived rules.
    """
    
    phrases = [
        'i ne znali',                            # and they did not know
        'i povelo že',                           # and it actually started leading
        'ne vymyla že',                          # and she did not actually wash thoroughly
        'ponesla',                               # she started carrying
        'vyvela že',                             # she actually led out
        'i unesli'                               # and they carried away
    ]

    def get_stress_syllable(phrase):
        """
        Applies the derived rules to a single phrase to find the stressed syllable.
        """
        vowels = "aeiouy"
        vowel_indices = [i for i, char in enumerate(phrase) if char in vowels]

        def get_syllable_from_char_index(char_index):
            """Finds the 1-based syllable number from a character's index in the phrase."""
            try:
                return vowel_indices.index(char_index) + 1
            except ValueError:
                # This should not happen if logic is correct
                return -1

        # Deconstruct the phrase
        words = phrase.split()
        has_i = 'i' in words
        has_ne = 'ne' in words
        has_že = 'že' in words
        
        # Isolate the verb form
        verb = ""
        for word in words:
            if word not in ['i', 'ne', 'že']:
                verb = word
                break

        # Rule 1: 'i ... -la' ending stress
        if has_i and verb.endswith('la'):
            stress_char_index = phrase.rfind('a')
            return get_syllable_from_char_index(stress_char_index)

        # Rule 2: 'ne' particle stress rules
        if has_ne:
            if verb.startswith('vy'):
                # Stress on 'vy-'
                stress_char_index = phrase.find('vy') + 1  # The 'y'
                return get_syllable_from_char_index(stress_char_index)
            else:
                # Stress on 'ne'
                stress_char_index = phrase.find('ne') + 1  # The 'e'
                return get_syllable_from_char_index(stress_char_index)

        # Rule 3: No 'ne', other rules
        if verb.startswith('vy'):
            # Stress on 'vy-'
            stress_char_index = phrase.find('vy') + 1  # The 'y'
            return get_syllable_from_char_index(stress_char_index)

        if has_že and verb.endswith('li'):
            # Stress on 'že'
            stress_char_index = phrase.rfind('že') + 1  # The 'e'
            return get_syllable_from_char_index(stress_char_index)

        # Rule 4: Default stem stress
        verb_start_in_phrase = phrase.find(verb)
        root_start_in_phrase = verb_start_in_phrase
        
        prefixes = ['po', 'u'] # Deduced from examples
        for p in prefixes:
            if verb.startswith(p):
                root_start_in_phrase += len(p)
                break
        
        # Find the first vowel of the root
        for idx in vowel_indices:
            if idx >= root_start_in_phrase:
                stress_char_index = idx
                return get_syllable_from_char_index(stress_char_index)
        
        return -1 # Should not be reached

    results = []
    for p in phrases:
        results.append(str(get_stress_syllable(p)))
    
    print("".join(results))

solve_old_russian_stress()