import re

def get_stress_position(phrase, stressed_component):
    """
    Calculates the 1-based syllable number of a stressed component in a phrase.
    A syllable is defined by a vowel (a, e, i, o, u, y).
    The position is determined by counting vowels up to and including the
    first vowel of the stressed component.
    If the component is 'ending', it counts all vowels in the phrase.
    """
    if stressed_component == 'ending':
        return len(re.findall('[aeiouy]', phrase))
    
    # Find the starting position of the stressed component
    try:
        start_index = phrase.find(stressed_component)
    except ValueError:
        return -1 # Should not happen

    # Count vowels in the part of the phrase leading up to the stressed component
    prefix_part = phrase[:start_index]
    vowel_count_before = len(re.findall('[aeiouy]', prefix_part))
    
    # The stressed syllable is the next one
    return vowel_count_before + 1

def solve_old_russian_stress():
    """
    Applies a derived set of linguistic rules to determine the stressed syllable
    in a list of Old Russian phrases and prints the results.
    """
    phrases_to_analyze = [
        {'text': 'i ne znali', 'root_type': 'attracting', 'root': 'zna'},
        {'text': 'i povelo že', 'root_type': 'non-attracting', 'root': 'vel'},
        {'text': 'ne vymyla že', 'root_type': 'attracting', 'root': 'myl'},
        {'text': 'ponesla', 'root_type': 'non-attracting', 'root': 'nes'},
        {'text': 'vyvela že', 'root_type': 'non-attracting', 'root': 'vel'},
        {'text': 'i unesli', 'root_type': 'non-attracting', 'root': 'nes'}
    ]

    final_digits = []
    
    for item in phrases_to_analyze:
        phrase = item['text']
        root_type = item['root_type']
        
        stressed_component = ''
        rule_explanation = ''

        # Apply the hierarchy of rules
        verb_part = phrase.replace('i ', '').replace('ne ', '').replace(' že', '')

        if verb_part.startswith('vy'):
            # Rule 1: 'vy-' prefix is always stressed.
            stressed_component = 'vy'
            rule_explanation = "Rule 1: The prefix 'vy-' is always stressed."
        elif root_type == 'non-attracting':
            # Rules for non-attracting roots
            if phrase.startswith('ne '):
                stressed_component = 'ne'
                rule_explanation = "Rule 2a: For a non-attracting root, an initial 'ne' is stressed."
            elif phrase.endswith(' že'):
                stressed_component = 'že'
                rule_explanation = "Rule 2b: For a non-attracting root, a final 'že' is stressed."
            else:
                stressed_component = 'ending'
                rule_explanation = "Rule 2c: For a non-attracting root, the ending is stressed by default."
        elif root_type == 'attracting':
            # Rule 3: Default for attracting roots is to stress the root.
            stressed_component = item['root']
            rule_explanation = "Rule 3: For a stress-attracting root, the root syllable is stressed."
        
        position = get_stress_position(phrase, stressed_component)
        final_digits.append(str(position))
        
        print(f"For '{phrase}':")
        print(f"  - {rule_explanation}")
        print(f"  - The stressed syllable is: {position}\n")

    print(f"The final sequence of stressed syllable numbers is: {''.join(final_digits)}")


solve_old_russian_stress()
<<<352314>>>