def solve_old_russian_stress():
    """
    This script determines the stressed syllable in six Old Russian phrases
    based on a set of linguistic rules derived from examples.
    """

    # The six phrases to be analyzed
    phrases_to_analyze = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    # Pre-analyzed information about the verbs and particles.
    # 'class': 'stem' or 'ending'
    # 'stress_on': for stem-stressed, the syllable that takes stress.
    # 'syllables': pre-calculated syllable count for each word.
    word_info = {
        # Particles
        'i': {'syllables': 1},
        'ne': {'syllables': 1},
        'že': {'syllables': 1},
        # Verbs
        'znali':   {'class': 'stem', 'stress_on': 'zna', 'syllables': 2},
        'povelo':  {'class': 'ending', 'syllables': 3},
        'vymyla':  {'class': 'stem', 'stress_on': 'vy', 'syllables': 3},
        'ponesla': {'class': 'ending', 'syllables': 3},
        'vyvela':  {'class': 'stem', 'stress_on': 'vy', 'syllables': 3},
        'unesli':  {'class': 'ending', 'syllables': 3},
    }

    # Morpheme-based syllables for stress location within a verb
    verb_syllables = {
        'znali': ['zna', 'li'],
        'vymyla': ['vy', 'my', 'la'],
        'vyvela': ['vy', 've', 'la'],
        'ponesla': ['po', 'nes', 'la'],
        'unesli': ['u', 'nes', 'li'],
    }

    final_results = []

    for phrase in phrases_to_analyze:
        words = phrase.split()
        verb = next(word for word in words if word in verb_syllables)
        info = word_info[verb]

        syllable_offset = 0
        stressed_syllable = 0

        if info['class'] == 'stem':
            # For stem-stressed verbs, find the position of the stressed stem syllable.
            for word in words:
                if word == verb:
                    # Find the index of the stressed syllable within the verb
                    stress_index_in_verb = verb_syllables[verb].index(info['stress_on'])
                    stressed_syllable = syllable_offset + 1 + stress_index_in_verb
                    break
                syllable_offset += word_info[word]['syllables']
        
        elif info['class'] == 'ending':
            # For ending-stressed verbs, apply rules based on particles and gender.
            total_syllables = sum(word_info[w]['syllables'] for w in words)
            
            if verb.endswith('a'):  # Feminine
                # Stress is on the final 'a' of the verb
                if 'že' in words:
                    stressed_syllable = total_syllables - 1
                else:
                    stressed_syllable = total_syllables
            else:  # Neuter or Plural
                if 'ne' in words:
                    # Stress is on 'ne'
                    stressed_syllable = word_info['i']['syllables'] + 1 if 'i' in words else 1
                elif 'že' in words:
                    # Stress is on 'že'
                    stressed_syllable = total_syllables
                else:
                    # Stress is on the verb's ending
                    stressed_syllable = total_syllables
        
        final_results.append(stressed_syllable)

    # Format the output as an equation, as requested.
    equation_parts = []
    for i, num in enumerate(final_results):
        power = len(final_results) - 1 - i
        equation_parts.append(f"{num}*10^{power}")
    
    final_number_str = "".join(map(str, final_results))
    final_number_int = int(final_number_str)

    print(f"The stressed syllables for the phrases are: {', '.join(map(str, final_results))}")
    print("The final number is formed by concatenating these digits.")
    print(f"Equation: {' + '.join(equation_parts)} = {final_number_int}")


solve_old_russian_stress()
<<<352314>>>