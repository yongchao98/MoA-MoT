import re

def count_vowels(s):
    """Counts the number of vowels in a string to determine its syllable count."""
    return len(re.findall(r'[aeiouy]', s))

def get_verb_parts(verb_word):
    """Parses a verb word to find its prefix and stem."""
    # Known stems and prefixes from the problem description
    stems = ['nes', 've', 'zna', 'my']
    prefixes = ['vy', 'po', 'u']

    found_stem = None
    for stem in stems:
        if stem in verb_word:
            found_stem = stem
            break

    found_prefix = None
    if found_stem:
        prefix_part = verb_word.split(found_stem)[0]
        if prefix_part in prefixes:
            found_prefix = prefix_part
            
    return found_prefix, found_stem

def find_stress_position(phrase):
    """
    Determines the stressed syllable in an Old Russian phrase based on a set of derived rules.
    """
    parts = phrase.split(' ')
    
    # 1. Identify particles and the main verb word
    has_i = parts[0] == 'i'
    has_ne = 'ne' in parts
    verb_word = ""
    for part in parts:
        if part not in ['i', 'ne', 'že']:
            verb_word = part
            break

    prefix, stem = get_verb_parts(verb_word)

    # 2. Determine the verb's stress class
    stress_class = ''
    if prefix == 'vy':
        stress_class = 'fixed_prefix'  # Rule 2a
    elif stem in ['nes', 've']:
        stress_class = 'movable'         # Class 1
    elif stem in ['zna', 'my']:
        stress_class = 'fixed_stem'      # Rule 2b

    # 3. Apply the rules for the determined class
    
    # For Movable Stress verbs (Class 1)
    if stress_class == 'movable':
        # Rule 1a: If phrase starts with 'ne', stress is on 'ne' (syllable 1)
        if has_ne and not has_i:
            return 1
        # Rule 1b: Otherwise, stress is on the last syllable of the entire phrase
        else:
            # The total number of vowels in the phrase gives the last syllable's position
            return count_vowels(phrase.replace(" ", ""))

    # For Fixed Stress verbs (Class 2)
    # Calculate the number of syllables before the verb word
    syllable_offset = 0
    if has_i: syllable_offset += 1
    if has_ne and parts[syllable_offset] == 'ne': syllable_offset += 1
    
    # Rule 2a: Stress on prefix 'vy-'
    if stress_class == 'fixed_prefix':
        # The prefix 'vy-' is the first syllable of the verb word
        return syllable_offset + 1
    
    # Rule 2b: Stress on the verb stem
    if stress_class == 'fixed_stem':
        # Count syllables in the prefix, if any (e.g., 'po-' in 'poznalo')
        prefix_syllables = count_vowels(prefix) if prefix else 0
        # Stress is after the prefix, on the stem
        return syllable_offset + prefix_syllables + 1
    
    return -1 # Should not be reached

# --- Main Execution ---

# The six phrases to be analyzed
phrases_to_analyze = [
    "i ne znali",
    "i povelo že",
    "ne vymyla že",
    "ponesla",
    "vyvela že",
    "i unesli"
]

print("Determining the stressed syllable for each phrase:")

results = []
for phrase in phrases_to_analyze:
    stress_position = find_stress_position(phrase)
    results.append(str(stress_position))
    # Output the result for each phrase
    print(f"'{phrase}': {stress_position}")

# Combine the results into a single string as requested
final_answer_string = "".join(results)
print("\nThe final answer is the sequence of the six numbers:")
print(final_answer_string)