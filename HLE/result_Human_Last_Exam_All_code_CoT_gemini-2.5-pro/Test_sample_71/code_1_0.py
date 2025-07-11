def get_stress_syllable(phrase: str) -> int:
    """
    Determines the stressed syllable number in an Old Russian phrase based on a set of derived rules.
    """
    vowels = "aeiouy"
    words = phrase.split()

    has_i = 'i' in words
    has_ne = 'ne' in words
    
    verb = ""
    for word in words:
        if word not in ['i', 'ne', 'že']:
            verb = word
            break

    # This variable will hold the character index of the stressed vowel in the phrase
    stressed_vowel_index = -1

    # Rule 1: Phrase contains 'i' and 'ne'. Stress is on the last syllable of the verb.
    if has_i and has_ne:
        verb_start_index = phrase.find(verb)
        # Find the index of the last vowel within the verb string
        last_vowel_pos_in_verb = -1
        for i, char in enumerate(verb):
            if char in vowels:
                last_vowel_pos_in_verb = i
        # Calculate the global index of that vowel in the full phrase
        stressed_vowel_index = verb_start_index + last_vowel_pos_in_verb

    # Rule 2: Phrase contains 'ne' (but not 'i').
    elif has_ne:
        # Rule 2a: If verb has 'vy-' prefix, stress is on the prefix.
        if verb.startswith('vy'):
            stressed_vowel_index = phrase.find('vy') + 1  # Index of 'y' in 'vy'
        # Rule 2b: Otherwise, stress is on 'ne'.
        else:
            stressed_vowel_index = phrase.find('ne') + 1  # Index of 'e' in 'ne'

    # Rule 3: Phrase does not contain 'ne'.
    else:
        is_unprefixed = not (verb.startswith('vy') or verb.startswith('po') or verb.startswith('u'))
        
        # Rule 3b: End-stress rule for specific unprefixed verbs.
        if (is_unprefixed and verb.endswith('li')) or \
           (is_unprefixed and verb.endswith('la') and has_i):
            # Find the last vowel in the entire phrase.
            for i in range(len(phrase) - 1, -1, -1):
                if phrase[i] in vowels:
                    stressed_vowel_index = i
                    break
        # Rule 3a: Stress on 'vy-' prefix.
        elif verb.startswith('vy'):
            stressed_vowel_index = phrase.find('vy') + 1 # Index of 'y'
        # Rule 3c: Default root-stress rule.
        else:
            verb_start_index = phrase.find(verb)
            vowel_indices_in_verb = [i for i, char in enumerate(verb) if char in vowels]
            if is_unprefixed:
                # First vowel of the verb (first syllable of the root)
                stressed_vowel_index = verb_start_index + vowel_indices_in_verb[0]
            else:  # Prefixed
                # Second vowel of the verb (first syllable of the root)
                stressed_vowel_index = verb_start_index + vowel_indices_in_verb[1]

    # Count the number of syllables (vowels) up to and including the stressed vowel.
    syllable_count = 0
    for i in range(stressed_vowel_index + 1):
        if phrase[i] in vowels:
            syllable_count += 1
            
    return syllable_count

# The list of phrases to analyze
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
print("Calculating the stressed syllable for each phrase:")
for phrase in phrases_to_solve:
    stress_position = get_stress_syllable(phrase)
    results.append(str(stress_position))
    print(f"'{phrase}': {stress_position}")

# Combine the results into a single string of digits
final_answer = "".join(results)
print("\nThe final combined answer is:")
print(final_answer)