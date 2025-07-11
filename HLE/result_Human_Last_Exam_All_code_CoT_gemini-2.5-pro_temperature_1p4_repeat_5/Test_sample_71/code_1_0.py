import re

def solve_old_russian_stress():
    """
    This script determines the stressed syllable in Old Russian phrases based on a derived set of rules.
    """

    # --- Step 1: Define the problem data ---
    # List of phrases to analyze
    phrases_to_analyze = [
        'i ne znali',
        'i povelo že',
        'ne vymyla že',
        'ponesla',
        'vyvela že',
        'i unesli'
    ]

    # Properties of verb roots derived from the examples
    mobile_roots = ['nes', 've']
    fixed_roots = ['zna', 'my']
    vowels = "aeiouy"
    
    # Store the final numerical results
    final_results = []

    print("Analyzing stressed syllables based on derived rules:\n")

    # --- Step 2: Process each phrase ---
    for phrase in phrases_to_analyze:
        
        # Prepare data for analysis
        words = phrase.split(' ')
        nospace_phrase = phrase.replace(" ", "")
        syllable_vowel_indices = [m.start() for m in re.finditer(f'[{vowels}]', nospace_phrase)]

        verb = ""
        for word in words:
            if word not in ['i', 'ne', 'že']:
                verb = word
                break
        
        stressed_syllable_num = 0

        # --- Step 3: Apply the rule hierarchy ---

        # Rule 1: 'vy-' prefix is always stressed.
        if verb.startswith('vy'):
            vy_vowel_index = nospace_phrase.find('vy') + 1
            stressed_syllable_num = syllable_vowel_indices.index(vy_vowel_index) + 1
        
        # Rule 2: Phrase starts with 'ne'.
        elif words[0] == 'ne':
            stressed_syllable_num = 1
            
        else:
            # Rule 3 & 4: Determine stress based on root type.
            root_class = None
            root_str = None
            for r in mobile_roots:
                if r in verb:
                    root_class = 'mobile'
                    break
            if not root_class:
                for r in fixed_roots:
                    if r in verb:
                        root_class = 'fixed'
                        root_str = r
                        break
            
            # Apply fixed-stress rule.
            if root_class == 'fixed':
                root_vowel = 'a' if root_str == 'zna' else 'y'
                # Find the vowel of the root in the phrase
                root_start_index = nospace_phrase.find(verb)
                verb_part = nospace_phrase[root_start_index:]
                root_vowel_index = root_start_index + verb_part.find(root_vowel)
                stressed_syllable_num = syllable_vowel_indices.index(root_vowel_index) + 1
            
            # Apply mobile-stress rule.
            elif root_class == 'mobile':
                # Stress falls on the last syllable of the phrase.
                stressed_syllable_num = len(syllable_vowel_indices)

        print(f"'{phrase}': The stressed syllable is number {stressed_syllable_num}")
        final_results.append(str(stressed_syllable_num))

    # --- Step 4: Output the final combined answer ---
    final_answer = "".join(final_results)
    print("\n------------------------------------")
    print(f"The six required digits in a row are: {final_answer}")
    print("------------------------------------")

# Execute the analysis
solve_old_russian_stress()
<<<352314>>>