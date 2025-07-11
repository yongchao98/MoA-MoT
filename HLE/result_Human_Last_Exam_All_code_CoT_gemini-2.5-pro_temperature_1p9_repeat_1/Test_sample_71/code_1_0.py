import re

def find_stressed_syllable(phrase):
    """
    This function determines the stressed syllable in an Old Russian phrase based on a set of derived rules.
    """

    # Helper function to find the indices of all vowels (a,e,i,o,u,y) in a string.
    # Each vowel corresponds to one syllable.
    def get_vowel_indices(p_str):
        # We can ignore 'že' as it has no vowel and doesn't affect syllable count.
        p_str = p_str.split(" že")[0]
        return [m.start() for m in re.finditer(r'[aeiouy]', p_str)]

    vowel_indices = get_vowel_indices(phrase)
    
    # Rule 1: The prefix 'vy-' is always stressed. This rule has the highest priority.
    # We check for verb stems that start with 'vy-'.
    if "vyvel" in phrase or "vynes" in phrase or "vymyl" in phrase:
        # Find the position of the vowel 'y' in 'vy-'.
        vy_vowel_pos = phrase.find("vy") + 1
        # Determine which syllable number this corresponds to.
        stress_syllable_num = vowel_indices.index(vy_vowel_pos) + 1
        return stress_syllable_num

    # Rule 2: Root-stressed verbs. The stress is always on the root syllable ('zna' or 'my').
    if "zna" in phrase or "my" in phrase:
        root_vowel_pos = -1
        if "zna" in phrase:
            # The stressed vowel is 'a' in 'zna'.
            root_vowel_pos = phrase.find("zna") + 2
        elif "my" in phrase:
            # The stressed vowel is 'y' in 'my'.
            root_vowel_pos = phrase.find("my") + 1
        stress_syllable_num = vowel_indices.index(root_vowel_pos) + 1
        return stress_syllable_num

    # Rule 3: Mobile-stress verbs (default case for roots like 'nes', 've').
    # Stress position depends on particles and grammatical endings.
    has_ne = phrase.startswith("ne ") or " ne " in phrase
    has_i = phrase.startswith("i ")
    clean_phrase = phrase.split(" že")[0]
    is_feminine = clean_phrase.endswith("a")

    if has_ne:
        if is_feminine:
            # For feminine forms with 'ne', stress moves to the final syllable.
            return len(vowel_indices)
        else:
            # For neuter or plural forms, stress moves to the particle 'ne'.
            ne_vowel_pos = phrase.find("ne") + 1
            stress_syllable_num = vowel_indices.index(ne_vowel_pos) + 1
            return stress_syllable_num

    if has_i:
        # With particle 'i' (and no 'ne'), stress moves to the final syllable.
        return len(vowel_indices)

    # Default case: No 'i' or 'ne' particle is present.
    if is_feminine:
        # For feminine forms, stress is on the first syllable of the verb form.
        # Since there are no preceding particles, this is syllable 1.
        return 1
    
    # This part should not be reached with the given test cases.
    return -1

# The list of phrases we need to analyze.
phrases_to_solve = [
    "i ne znali",                              # Root 'zna' -> root-stressed
    "i povelo že",                             # Root 've' (mobile) + 'i' -> ending stress
    "ne vymyla že",                            # Prefix 'vy-' -> vy-stress
    "ponesla",                                 # Root 'nes' (mobile), feminine, no particles -> 1st syllable stress
    "vyvela že",                               # Prefix 'vy-' -> vy-stress
    "i unesli"                                 # Root 'nes' (mobile) + 'i' -> ending stress
]

# Calculate the stressed syllable for each phrase.
r1 = find_stressed_syllable(phrases_to_solve[0])
r2 = find_stressed_syllable(phrases_to_solve[1])
r3 = find_stressed_syllable(phrases_to_solve[2])
r4 = find_stressed_syllable(phrases_to_solve[3])
r5 = find_stressed_syllable(phrases_to_solve[4])
r6 = find_stressed_syllable(phrases_to_solve[5])

# Print the result for each phrase and the final combined number sequence.
print(f"The stressed syllable for 'i ne znali' is: {r1}")
print(f"The stressed syllable for 'i povelo že' is: {r2}")
print(f"The stressed syllable for 'ne vymyla že' is: {r3}")
print(f"The stressed syllable for 'ponesla' is: {r4}")
print(f"The stressed syllable for 'vyvela že' is: {r5}")
print(f"The stressed syllable for 'i unesli' is: {r6}")
print("\nFinal Answer:")
print(f"{r1}{r2}{r3}{r4}{r5}{r6}")
