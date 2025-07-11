import re

def get_stressed_syllable(phrase):
    """
    Analyzes an Old Russian phrase to find the stressed syllable based on a derived set of rules.
    """
    
    # ---- 1. Parse the phrase into its components ----
    
    parts = {'i': False, 'ne': False, 'ze': False, 'prefix': None, 'root': None, 'ending': None}
    original_phrase = phrase
    
    words = phrase.split(' ')

    # Particles and conjunctions
    if words[0] == 'i':
        parts['i'] = True
        words.pop(0)
    if words[0] == 'ne':
        parts['ne'] = True
        words.pop(0)
    if words[-1] == 'že':
        parts['ze'] = True
        words.pop(-1)
    
    verb = words[0]

    # Verb ending
    if verb.endswith('la'):
        parts['ending'] = 'la'
        verb_base = verb[:-2]
    elif verb.endswith('lo'):
        parts['ending'] = 'lo'
        verb_base = verb[:-2]
    elif verb.endswith('li'):
        parts['ending'] = 'li'
        verb_base = verb[:-2]
    else:
        # Should not happen with the given data
        verb_base = verb

    # Verb prefix and root
    prefixes = ['vy', 'po', 'u']
    for p in prefixes:
        if verb_base.startswith(p):
            parts['prefix'] = p
            parts['root'] = verb_base[len(p):]
            break
    
    if not parts['prefix']:
        parts['root'] = verb_base
        
    # ---- 2. Apply stress rules to find the stressed component ----
    
    root_stressed_roots = ['zna', 'my']
    stress_location = ''
    
    # Rule 1: 'vy-' prefix is always stressed.
    if parts['prefix'] == 'vy':
        stress_location = 'prefix'
    # Rule 2: 'ne' is present.
    elif parts['ne']:
        if parts['ending'] == 'la':
            stress_location = 'ending'
        else:
            stress_location = 'ne'
    # Rule 4 (Special Case): for phrases like 'nesli že'
    elif not parts['prefix'] and parts['root'] == 'nes' and parts['ending'] == 'li' and parts['ze']:
        stress_location = 'ze'
    # Rule 3: Default root-based or end-based stress.
    else:
        if parts['root'] in root_stressed_roots:
            stress_location = 'root'
        else: # Assumed end-stressed roots ('nes', 've')
            stress_location = 'ending'
            
    # ---- 3. Find the syllable number of the stressed component ----
    
    vowels = "aeiouy"
    stressed_substring = ""
    
    # Determine which part of the string to search for the vowel
    if stress_location == 'prefix':
        stressed_substring = parts['prefix']
        search_area = verb
    elif stress_location == 'ne':
        stressed_substring = 'ne'
        search_area = 'ne'
    elif stress_location == 'ze':
        stressed_substring = 'že'
        search_area = 'že'
    elif stress_location == 'root':
        stressed_substring = parts['root']
        search_area = verb
    elif stress_location == 'ending':
        stressed_substring = parts['ending']
        search_area = verb
    
    # Find the position of the stressed vowel within the whole phrase
    stress_vowel_abs_pos = -1
    # Find start of the component containing the stress
    component_start = original_phrase.find(search_area)
    # Find start of the specific stressed part within that component
    substring_start = search_area.find(stressed_substring)
    
    # Find the first vowel in the stressed substring
    for i, char in enumerate(stressed_substring):
        if char in vowels:
            stress_vowel_abs_pos = component_start + substring_start + i
            break
            
    # Count vowels from the beginning of the phrase to find syllable number
    syllable_count = 0
    for i, char in enumerate(original_phrase):
        if char in vowels:
            syllable_count += 1
            if i == stress_vowel_abs_pos:
                return syllable_count
    
    return -1 # Should not be reached

if __name__ == '__main__':
    phrases = [
        'i ne znali',
        'i povelo že',
        'ne vymyla že',
        'ponesla',
        'vyvela že',
        'i unesli'
    ]
    
    results = []
    print("Calculating stressed syllable for each phrase:")
    for phrase in phrases:
        result = get_stressed_syllable(phrase)
        print(f"{phrase} -> {result}")
        results.append(str(result))
    
    final_answer = "".join(results)
    print(f"\nThe final sequence of digits is {final_answer}.")
    print(f"<<<{final_answer}>>>")
