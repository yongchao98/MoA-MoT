import re

def get_stress_syllable(phrase):
    """
    Calculates the stressed syllable number for an Old Russian phrase based on derived rules.
    """
    phrase = phrase.strip()
    
    # 1. Find all vowels and their indices to count syllables.
    vowels = "aeiouy"
    vowel_indices = [i for i, char in enumerate(phrase) if char in vowels]
    
    # 2. Identify the verb root to determine its stress group.
    group1_roots = ['nes', 'nos', 'vel', 'ved']
    group2_roots = ['zna', 'my']
    
    root = ""
    root_group = 0
    root_start_index = -1

    for r in group1_roots:
        if r in phrase:
            root = r
            root_group = 1
            root_start_index = phrase.find(r)
            break
    if not root:
        for r in group2_roots:
            if r in phrase:
                root = r
                root_group = 2
                root_start_index = phrase.find(r)
                break

    # 3. Apply stress rules based on the verb group.

    # Group 2: Fixed stress on the root syllable.
    if root_group == 2:
        # Find the index of the vowel within the root substring.
        root_vowel_index_in_phrase = -1
        for i, char in enumerate(root):
            if char in vowels:
                root_vowel_index_in_phrase = root_start_index + i
                break
        # Find which syllable number this vowel corresponds to.
        stress_syllable_num = vowel_indices.index(root_vowel_index_in_phrase) + 1
        return stress_syllable_num
        
    # Group 1: Mobile stress with hierarchical rules.
    elif root_group == 1:
        # Rule a: `vy-` prefix is always stressed.
        # It can be at the start of the phrase or after a particle like 'ne'.
        vy_match = re.search(r'(^|\s)vy', phrase)
        if vy_match:
            # The stressed vowel is 'y' in 'vy'.
            vy_vowel_index = vy_match.start(0) + vy_match.group(0).find('y')
            stress_syllable_num = vowel_indices.index(vy_vowel_index) + 1
            return stress_syllable_num

        # Rule b: `i` particle shifts stress to the end. This rule takes precedence over 'ne'.
        if phrase.startswith('i '):
            # Stress falls on the last syllable of the phrase.
            return len(vowel_indices)

        # Rule c: `ne` particle takes the stress (if not overridden by `vy-` or `i`).
        if phrase.startswith('ne '):
            # Stress falls on the 'e' of 'ne'.
            ne_vowel_index = phrase.find('e')
            stress_syllable_num = vowel_indices.index(ne_vowel_index) + 1
            return stress_syllable_num
            
        # Rule d: Default for Group 1 is stress on the last syllable.
        return len(vowel_indices)

    return 0 # Should not be reached

# The six phrases to be analyzed.
phrases_to_analyze = [
    'i ne znali',                            # Expected: 3
    'i povelo že',                           # Expected: 5
    'ne vymyla že',                          # Expected: 3
    'ponesla',                               # Expected: 3
    'vyvela že',                             # Expected: 1
    'i unesli'                               # Expected: 4
]

# Calculate the stressed syllable for each phrase and build the final answer string.
result_digits = []
for p in phrases_to_analyze:
    stress_pos = get_stress_syllable(p)
    result_digits.append(str(stress_pos))

final_answer = "".join(result_digits)
print(final_answer)