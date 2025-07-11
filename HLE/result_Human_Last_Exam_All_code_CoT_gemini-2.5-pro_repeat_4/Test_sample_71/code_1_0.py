def solve_stress():
    """
    This function analyzes Old Russian phrases to determine the stressed syllable.
    It implements a rule-based system based on verb classification and particles.
    """

    def get_syllable_number(phrase, index):
        """
        Calculates the syllable number for a given character index in a phrase.
        Each vowel is considered to mark a syllable.
        """
        vowels = "aeiouy"
        syllable_count = 0
        # Iterate through the phrase up to the stressed character
        # and count the vowels to find the syllable number.
        for i in range(index + 1):
            if phrase[i] in vowels:
                syllable_count += 1
        return syllable_count

    def calculate_stress_position(phrase):
        """
        Applies a set of linguistic rules to find the stressed syllable in a phrase.
        """
        words = phrase.split()
        
        # Identify particles and the main verb
        has_i = "i" in words
        has_ne = "ne" in words
        has_ze = "že" in words
        # The verb is whatever is left after removing particles
        verb_word = [w for w in words if w not in ["i", "ne", "že"]][0]

        # 1. Identify verb class based on its root
        root_class = None
        if "nes" in verb_word or "ve" in verb_word:
            root_class = "A"  # Mobile stress
        elif "zna" in verb_word or "my" in verb_word:
            root_class = "B"  # Fixed root stress

        stress_index = -1

        # 2. Apply rules based on verb class
        
        # Rule for Class B verbs: Stress is always on the root
        if root_class == "B":
            root_substring = "zna" if "zna" in verb_word else "my"
            stress_index = phrase.find(root_substring)

        # Rules for Class A verbs (in order of priority)
        elif root_class == "A":
            # Priority 1: 'ne' is present AND 'že' is absent -> stress on 'ne'
            if has_ne and not has_ze:
                stress_index = phrase.find("ne")
            # Priority 2: 'i' is present -> stress on the verb's ending
            elif has_i:
                verb_start_index = phrase.find(verb_word)
                verb_end_index = verb_start_index + len(verb_word)
                # Find the index of the last vowel in the verb word
                for i in range(verb_start_index, verb_end_index):
                    if phrase[i] in "aeiouy":
                        stress_index = i
            # Priority 3: 'že' is present (and 'i' is absent)
            elif has_ze:
                # If prefix is 'vy-', stress is on the prefix
                if verb_word.startswith("vy"):
                    stress_index = phrase.find("vy")
                # Otherwise (unprefixed), stress is on the ending
                else:
                    verb_start_index = phrase.find(verb_word)
                    verb_end_index = verb_start_index + len(verb_word)
                    for i in range(verb_start_index, verb_end_index):
                        if phrase[i] in "aeiouy":
                            stress_index = i
            # Default: No special particles -> stress is on the root
            else:
                root_substring = "nes" if "nes" in verb_word else "ve"
                stress_index = phrase.find(root_substring)
        
        return get_syllable_number(phrase, stress_index)

    phrases_to_analyze = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli",
    ]

    final_result = ""
    for p in phrases_to_analyze:
        final_result += str(calculate_stress_position(p))

    print(final_result)

solve_stress()