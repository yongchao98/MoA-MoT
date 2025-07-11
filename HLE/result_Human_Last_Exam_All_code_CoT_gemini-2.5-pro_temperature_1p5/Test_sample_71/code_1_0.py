def solve_old_russian_stress():
    """
    This script deduces and applies stress rules for Old Russian phrases.
    """
    phrases = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    # --- Data Definitions ---
    prefixes_set = {'vy', 'po', 'u'}
    roots_map = {
        'nes': 'E', 've': 'E',  # Ending-stressed group
        'zna': 'R', 'my': 'R'   # Root-stressed group
    }
    endings_set = {'la', 'lo', 'li'}
    vowels = "aeiouy"

    def parse_phrase(phrase):
        """Parses a phrase into its grammatical components."""
        components = {
            'pre_particles': [], 'post_particles': [],
            'prefix': '', 'root': '', 'ending': '', 'root_group': '',
            'verb_word': '', 'full_phrase': phrase
        }
        words = phrase.split()
        
        # Identify particles and verb
        verb_word_found = ""
        for word in words:
            if word in ['i', 'ne']:
                components['pre_particles'].append(word)
            elif word == 'že':
                components['post_particles'].append(word)
            else:
                verb_word_found = word
                components['verb_word'] = word

        # Parse the verb word
        temp_verb = verb_word_found
        for p in prefixes_set:
            if temp_verb.startswith(p):
                components['prefix'] = p
                temp_verb = temp_verb[len(p):]
                break
        for r, r_type in roots_map.items():
            if temp_verb.startswith(r):
                components['root'] = r
                components['root_group'] = r_type
                temp_verb = temp_verb[len(r):]
                break
        components['ending'] = temp_verb
        
        return components

    def get_stressed_part_and_rule(c):
        """Applies the deduced stress rules."""
        # Rule 1: 'vy-' is always stressed
        if c['prefix'] == 'vy':
            return c['prefix'], "Rule 1: 'vy-' prefix is always stressed."

        # Rule 2: 'ne' + root 'nes-'
        if 'ne' in c['pre_particles'] and c['root'] == 'nes':
            if c['ending'] in ['lo', 'li']:
                return 'ne', "Rule 2a: 'ne' with root 'nes-' and ending '-lo/-li' stresses 'ne'."
            if c['ending'] == 'la':
                return c['ending'], "Rule 2b: 'ne' with root 'nes-' and ending '-la' stresses the ending."
        
        # Rule 3: Unprefixed Ending-stressed group (no 'ne')
        if not c['prefix'] and c['root_group'] == 'E' and 'ne' not in c['pre_particles']:
             rule_desc = "Rule 3: Unprefixed ending-stressed verb stresses the last syllable of the phrase."
             if c['post_particles']:
                 return c['post_particles'][-1], rule_desc
             else:
                 return c['ending'], rule_desc

        # Rule 4: Default is root stress
        return c['root'], "Rule 4: Default stress is on the root."

    def count_stress_position(full_phrase, components, stressed_part):
        """Counts the syllable number for the stressed part."""
        phrase_for_count = full_phrase.replace('že', 'ze')
        
        # Determine the word where the stress falls
        target_word = ""
        if stressed_part in ['i', 'ne', 'že', 'ze']:
            target_word = stressed_part
        else:
            target_word = components['verb_word']

        # Find the starting position of the stressed morpheme
        word_start_idx = phrase_for_count.find(target_word)
        morpheme_start_idx = target_word.find(stressed_part)
        # For endings, find the last occurrence
        if morpheme_start_idx == -1 and stressed_part == components['ending']:
            morpheme_start_idx = target_word.rfind(stressed_part)
            
        abs_stress_point = word_start_idx + morpheme_start_idx

        # Count vowels up to the start of the stressed morpheme
        syllable_num = 0
        for i in range(abs_stress_point):
            if phrase_for_count[i] in vowels:
                syllable_num += 1
        
        return syllable_num + 1

    final_digits = []
    print("Determining stress for each phrase:")
    for phrase in phrases:
        comps = parse_phrase(phrase)
        stressed_part, rule = get_stressed_part_and_rule(comps)
        position = count_stress_position(phrase, comps, stressed_part)
        final_digits.append(str(position))
        
        print(f"- '{phrase}': Stress on '{stressed_part}' (Reason: {rule}) -> Syllable number: {position}")

    final_answer = "".join(final_digits)
    print("\nThe final answer is the sequence of these six digits.")
    print(final_answer)
    return final_answer

final_answer_string = solve_old_russian_stress()
print(f"<<<{final_answer_string}>>>")