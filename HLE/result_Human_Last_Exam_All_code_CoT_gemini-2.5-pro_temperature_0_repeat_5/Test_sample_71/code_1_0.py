def find_stress_in_old_russian():
    """
    This script determines the stressed syllable in Old Russian phrases based on a set of derived rules.
    """
    phrases = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    results = []
    print("Determining the stressed syllable for each phrase:")
    print("-" * 50)

    for phrase in phrases:
        # 1. Parse the phrase into components and an ordered list of syllables (morphemes)
        parts = phrase.split(' ')
        
        # Components dictionary for rule checking
        components = {
            'has_i': 'i' in parts,
            'has_ne': 'ne' in parts,
            'has_zhe': 'že' in parts,
            'prefix': None,
            'root': None,
            'ending': None,
            'verb': None
        }
        
        # Ordered list of morphemes for syllable counting
        syllable_parts = []
        
        # Known morphemes
        PREFIXES = ['vy', 'po', 'u']
        ENDINGS = ['la', 'lo', 'li']
        PARTICLES = ['i', 'ne', 'že']

        for part in parts:
            if part in PARTICLES:
                syllable_parts.append({'type': 'particle', 'text': part})
            else: # It's a verb
                components['verb'] = part
                verb_to_parse = part
                
                # Find prefix
                for p in PREFIXES:
                    if verb_to_parse.startswith(p):
                        components['prefix'] = p
                        syllable_parts.append({'type': 'prefix', 'text': p})
                        verb_to_parse = verb_to_parse[len(p):]
                        break
                
                # The rest is root + ending. Find ending from the end.
                temp_root = verb_to_parse
                for e in ENDINGS:
                    if verb_to_parse.endswith(e):
                        components['ending'] = e
                        temp_root = verb_to_parse[:-len(e)]
                        break
                
                components['root'] = temp_root
                syllable_parts.append({'type': 'root', 'text': temp_root})

                if components['ending']:
                    syllable_parts.append({'type': 'ending', 'text': components['ending']})

        # 2. Apply stress rules
        stressed_part_type = None
        stressed_part_text = None
        rule_applied = ""

        # Rule 1: 'vy-' prefix is always stressed
        if components['prefix'] == 'vy':
            stressed_part_type = 'prefix'
            stressed_part_text = 'vy'
            rule_applied = "Rule 1: Prefix 'vy-' is stressed."
        # Rule 2: 'ne' is stressed (if 'i' is not present)
        elif components['has_ne'] and not components['has_i']:
            stressed_part_type = 'particle'
            stressed_part_text = 'ne'
            rule_applied = "Rule 2: Particle 'ne' is stressed."
        # Rule 3: 'i' + '-la' ending -> ending is stressed
        elif components['has_i'] and components['ending'] == 'la':
            stressed_part_type = 'ending'
            stressed_part_text = 'la'
            rule_applied = "Rule 3: With 'i' and ending '-la', the ending is stressed."
        # Rule 4: 'že' + prefixless verb + '-li' ending -> 'že' is stressed
        elif components['has_zhe'] and components['prefix'] is None and components['ending'] == 'li':
            stressed_part_type = 'particle'
            stressed_part_text = 'že'
            rule_applied = "Rule 4: With 'že', no prefix, and ending '-li', 'že' is stressed."
        # Rule 5: Default to root stress
        else:
            stressed_part_type = 'root'
            stressed_part_text = components['root']
            rule_applied = "Rule 5: Default stress on the root."

        # 3. Find the syllable number
        stress_pos = -1
        for i, part in enumerate(syllable_parts):
            if part['type'] == stressed_part_type:
                # For root, type is enough. For others, check text to be sure.
                if stressed_part_type == 'root' or part['text'] == stressed_part_text:
                    stress_pos = i + 1
                    break
        
        results.append(str(stress_pos))
        print(f"Phrase: '{phrase}'")
        print(f"Analysis: {rule_applied}")
        print(f"Stressed Syllable: {stress_pos} (morpheme: '{stressed_part_text}')")
        print("-" * 50)

    final_answer = "".join(results)
    print(f"The final sequence of stressed syllables is: {final_answer}")

if __name__ == '__main__':
    find_stress_in_old_russian()