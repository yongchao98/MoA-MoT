def solve_old_russian_stress():
    """
    Determines the stressed syllable in Old Russian phrases by applying a set of derived rules.
    The rules are based on analyzing patterns in the provided examples.
    """

    # The list of phrases for which we need to find the stress.
    phrases_to_solve = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    final_results = []

    # Helper function to count the number of syllables up to a specific character index.
    # In this Romanization, each vowel corresponds to one syllable.
    def count_syllables_to_index(phrase, index):
        vowels = "aeiouy"
        syllable_count = 0
        for i, char in enumerate(phrase):
            if i > index:
                break
            if char in vowels:
                syllable_count += 1
        return syllable_count

    # Helper function to find the index of the Nth vowel within a given word.
    def find_nth_vowel_index_in_word(word, n):
        vowels = "aeiouy"
        vowel_indices = [i for i, char in enumerate(word) if char in vowels]
        if n > 0 and len(vowel_indices) >= n:
            return vowel_indices[n - 1]
        return -1

    for phrase in phrases_to_solve:
        # --- Rule Application ---

        # Rule 1: The prefix 'vy-' is always stressed. This rule has the highest priority.
        if 'vy' in phrase:
            # The stress falls on the 'y' of 'vy'.
            vy_index = phrase.find('vy')
            final_results.append(count_syllables_to_index(phrase, vy_index + 1))
            continue

        # Rule 2: The negative particle 'ne' is usually stressed.
        if 'ne ' in phrase:
            # Exception: If the verb ends in '-a', the ending is stressed instead.
            if phrase.endswith('a'):
                stress_index = phrase.rfind('a')
                final_results.append(count_syllables_to_index(phrase, stress_index))
            # Otherwise, 'ne' itself is stressed.
            else:
                stress_index = phrase.find('ne ') + 1  # Stress on the 'e' of 'ne'
                final_results.append(count_syllables_to_index(phrase, stress_index))
            continue

        # Rule 3: For all other cases (no 'vy-' and no 'ne' particle).
        # We need to identify the main verb and its properties.
        words = phrase.split(' ')
        verb_word = words[-2] if words[-1] == 'že' else words[-1]
        verb_start_index = phrase.find(verb_word)
        prefixes = ['po', 'u']
        has_prefix = any(verb_word.startswith(p) for p in prefixes)

        # Rule 3a: Verb ends in '-a'.
        if verb_word.endswith('a'):
            # Stress is on the ending '-a' if a prefix or the particle 'i' is present.
            if has_prefix or phrase.startswith('i '):
                stress_index = phrase.rfind('a')
                final_results.append(count_syllables_to_index(phrase, stress_index))
            # Otherwise, stress is on the root (the first syllable of the verb).
            else:
                verb_vowel_index = find_nth_vowel_index_in_word(verb_word, 1)
                final_results.append(count_syllables_to_index(phrase, verb_start_index + verb_vowel_index))

        # Rule 3b: Verb ends in '-o'.
        elif verb_word.endswith('o'):
            # Stress is on the root. The root is the 2nd syllable if there's a prefix,
            # and the 1st syllable if there is no prefix.
            syllable_to_stress = 2 if has_prefix else 1
            verb_vowel_index = find_nth_vowel_index_in_word(verb_word, syllable_to_stress)
            final_results.append(count_syllables_to_index(phrase, verb_start_index + verb_vowel_index))

        # Rule 3c: Verb ends in '-i'.
        elif verb_word.endswith('i'):
            # If a prefix is present, stress is on the root (2nd syllable of the verb).
            if has_prefix:
                verb_vowel_index = find_nth_vowel_index_in_word(verb_word, 2)
                final_results.append(count_syllables_to_index(phrase, verb_start_index + verb_vowel_index))
            # Otherwise, stress is on the ending '-i'.
            else:
                stress_index = phrase.rfind('i')
                final_results.append(count_syllables_to_index(phrase, stress_index))

    # Print the final result, which is the sequence of six digits.
    print("".join(map(str, final_results)))

solve_old_russian_stress()