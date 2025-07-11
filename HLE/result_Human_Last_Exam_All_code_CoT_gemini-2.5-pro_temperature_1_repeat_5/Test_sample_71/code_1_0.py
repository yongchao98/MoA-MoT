import re

def find_vowel_indices(text):
    """Finds the character indices of all vowels in a string."""
    return [m.start() for m in re.finditer(r'[aeiouy]', text)]

def get_verb_and_prefix(phrase):
    """Identifies the verb and its prefix in a phrase."""
    words = phrase.split()
    verb_word = ""
    for w in words:
        if 'l' in w:
            verb_word = w
            break
            
    prefix = None
    if verb_word.startswith('vy'): prefix = 'vy'
    elif verb_word.startswith('po'): prefix = 'po'
    elif verb_word.startswith('u'): prefix = 'u'
    
    return verb_word, prefix

def get_stressed_syllable(phrase):
    """
    Determines the stressed syllable number based on a set of hierarchical rules.
    Each vowel is considered a syllable.
    """
    has_ze = 'že' in phrase.split()
    has_ne = 'ne' in phrase.split()
    verb, prefix = get_verb_and_prefix(phrase)
    
    # Get the character position of every vowel in the phrase
    vowel_indices = find_vowel_indices(phrase)
    
    # --- Rule 1: The prefix 'vy-' is always stressed. ---
    if prefix == 'vy':
        # Find the position of 'y' in 'vy'
        vy_pos = phrase.find('vy')
        stress_char_index = vy_pos + 1
        # Return the syllable number (1-based index)
        return vowel_indices.index(stress_char_index) + 1
        
    # --- Rule 2: Feminine ending '-a' is stressed if 'že' is NOT present. ---
    if verb.endswith('a') and not has_ze:
        # Find the position of the last 'a' in the verb
        stress_char_index = phrase.rfind(verb) + len(verb) - 1
        return vowel_indices.index(stress_char_index) + 1
        
    # --- Rule 3: The particle 'ne' takes the stress. ---
    # This rule is overridden by rules 1 and 2.
    if has_ne:
        # Find the position of 'e' in 'ne'
        ne_pos = phrase.find('ne ')
        stress_char_index = ne_pos + 1
        return vowel_indices.index(stress_char_index) + 1
        
    # --- Rule 4: 'že' is stressed with unprefixed, plural verbs. ---
    if prefix is None and verb.endswith('i') and has_ze:
        # Find the position of 'e' in 'že'
        ze_pos = phrase.find('že')
        stress_char_index = ze_pos + 1
        return vowel_indices.index(stress_char_index) + 1
        
    # --- Rule 5: Default stress is on the verb's root. ---
    # The root is the first vowel of an unprefixed verb, or the second vowel of a prefixed verb.
    verb_vowel_indices = find_vowel_indices(verb)
    verb_pos_in_phrase = phrase.find(verb)
    
    if prefix:
        # Root is the second vowel of the verb
        root_vowel_pos_in_verb = verb_vowel_indices[1]
    else:
        # Root is the first vowel of the verb
        root_vowel_pos_in_verb = verb_vowel_indices[0]
        
    stress_char_index = verb_pos_in_phrase + root_vowel_pos_in_verb
    return vowel_indices.index(stress_char_index) + 1

# List of phrases to analyze
phrases_to_solve = [
    "i ne znali",
    "i povelo že",
    "ne vymyla že",
    "ponesla",
    "vyvela že",
    "i unesli"
]

# Calculate and store the results
results = []
print("Determining the stressed syllable for each phrase:")
for p in phrases_to_solve:
    stress_position = get_stressed_syllable(p)
    results.append(str(stress_position))
    print(f"'{p}' -> {stress_position}")

# Print the final answer in the required format
final_answer = "".join(results)
print("\nThe final sequence of digits is:")
print(final_answer)
