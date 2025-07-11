def solve_old_russian_stress():
    """
    Solves the Old Russian stress puzzle by applying a set of derived rules.
    """
    
    phrases = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    final_result_digits = []

    for phrase in phrases:
        # Step 1: Parse the phrase into its morphemes (which act as syllables)
        words = phrase.split(' ')
        syllables = []
        
        # Information for rule application
        verb_word = ""
        prefix = None
        root = None
        ending = None
        has_i = 'i' in words
        has_ne = 'ne' in words
        has_že = 'že' in words
        is_ne_initial = phrase.startswith('ne ')

        # Add initial particles to syllables list
        if has_i:
            syllables.append('i')
        if has_ne:
            syllables.append('ne')
            
        # Find the verb word
        for w in words:
            if w not in ['i', 'ne', 'že']:
                verb_word = w
                break
        
        temp_verb = verb_word
        
        # Parse prefix
        if temp_verb.startswith('vy'):
            prefix = 'vy'
            temp_verb = temp_verb[2:]
        elif temp_verb.startswith('po'):
            prefix = 'po'
            temp_verb = temp_verb[2:]
        elif temp_verb.startswith('u'):
            prefix = 'u'
            temp_verb = temp_verb[1:]
        
        if prefix:
            syllables.append(prefix)
            
        # Parse root and ending
        # The four known roots are: 'nes', 'zna', 'my', 've'
        if temp_verb.startswith('nes'):
            root = 'nes'
            ending = temp_verb[3:]
        elif temp_verb.startswith('zna'):
            root = 'zna'
            ending = temp_verb[3:]
        elif temp_verb.startswith('my'):
            root = 'my'
            ending = temp_verb[2:]
        elif temp_verb.startswith('ve'):
            root = 've'
            ending = temp_verb[2:]

        syllables.append(root)
        syllables.append(ending)
        
        # Add final particle
        if has_že:
            syllables.append('že')
            
        # Step 2: Apply the rules in order of priority to find the stressed syllable
        stress_pos = -1

        # Rule 1: Prefix `vy-`
        if prefix == 'vy':
            stress_pos = syllables.index('vy') + 1
        # Rule 2: Initial `ne`
        elif is_ne_initial:
            stress_pos = syllables.index('ne') + 1
        # Rule 3: `i ne ...`
        elif has_i and has_ne:
            stress_pos = syllables.index(ending) + 1
        # Rule 4: Unprefixed plural + `že`
        elif prefix is None and ending in ['li', 'i'] and has_že:
            stress_pos = syllables.index('že') + 1
        # Rule 5: Unprefixed fem. root `ve`
        elif prefix is None and root == 've' and ending in ['la', 'a']:
            stress_pos = syllables.index(ending) + 1
        # Rule 6: Default to root stress
        else:
            stress_pos = syllables.index(root) + 1
            
        final_result_digits.append(str(stress_pos))

    # Print the final result as a single string of digits
    print("".join(final_result_digits))

solve_old_russian_stress()