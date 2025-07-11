def solve():
    """
    Analyzes Old Russian phrases to determine the stressed syllable based on derived rules.
    """

    def get_stress_position(phrase, target_substring):
        """Finds the syllable number of the first vowel in the target substring."""
        vowels = "aeiouy"
        pos = phrase.find(target_substring)
        if pos == -1: return -1  # Should not happen
        
        syllable_count = 0
        # Count vowels up to the target substring's position
        for i in range(pos):
            if phrase[i] in vowels:
                syllable_count += 1
        
        # The target's first vowel constitutes the next syllable
        return syllable_count + 1

    def get_last_syllable_position(phrase, subphrase):
        """Finds the syllable number corresponding to the last vowel in the subphrase."""
        vowels = "aeiouy"
        # Find the portion of the main phrase that constitutes the subphrase
        pos = phrase.find(subphrase)
        end_pos = pos + len(subphrase)
        
        syllable_count = 0
        # Count all vowels up to the end of the subphrase
        for i in range(end_pos):
            if phrase[i] in vowels:
                syllable_count += 1
        return syllable_count

    def get_stressed_syllable(phrase):
        """Applies the derived rules to find the stressed syllable."""
        words = phrase.split(' ')
        verb = ""
        verb_index = -1
        for i, word in enumerate(words):
            if word not in ["i", "ne", "že"]:
                verb = word
                verb_index = i
                break
        
        # Rule 1: Verb is prefixed with 'vy-'
        if verb.startswith('vy'):
            return get_stress_position(phrase, 'vy')

        # Rule 2: Root is 'zna' or 'myl'
        if 'zna' in verb or 'myl' in verb:
            root = 'zna' if 'zna' in verb else 'myl'
            return get_stress_position(phrase, root)
            
        # Rule 3: Root is from the mobile stress group ('nes', 'vel')
        has_ze = "že" in words
        has_ne = "ne" in words
        has_i = "i" in words
        
        if has_ze:
            # Stress is on 'že'
            return get_stress_position(phrase, "že")
        
        if has_ne:
            if has_i:
                # With 'i ne', stress is on the verb's final syllable
                verb_phrase_segment = " ".join(words[:verb_index + 1])
                return get_last_syllable_position(phrase, verb_phrase_segment)
            else:
                # With 'ne' alone, stress is on 'ne'
                return 1
                
        # With 'i' alone or no particles, stress is on the verb's final syllable
        verb_phrase_segment = " ".join(words[:verb_index + 1])
        return get_last_syllable_position(phrase, verb_phrase_segment)

    # --- Main execution ---
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

    # Phrase 1
    p1 = phrases[0]
    r1 = get_stressed_syllable(p1)
    results.append(r1)
    print(f"Phrase: '{p1}'")
    print(f"Rule: Root is 'zna' (Group 2), so stress is on the root.")
    print(f"Result: The stressed syllable is number {r1}\n")

    # Phrase 2
    p2 = phrases[1]
    r2 = get_stressed_syllable(p2)
    results.append(r2)
    print(f"Phrase: '{p2}'")
    print(f"Rule: Mobile stress verb with particle 'že' (Group 3), so stress is on 'že'.")
    print(f"Result: The stressed syllable is number {r2}\n")

    # Phrase 3
    p3 = phrases[2]
    r3 = get_stressed_syllable(p3)
    results.append(r3)
    print(f"Phrase: '{p3}'")
    print(f"Rule: Verb has prefix 'vy-' (Group 1), so stress is on 'vy'.")
    print(f"Result: The stressed syllable is number {r3}\n")

    # Phrase 4
    p4 = phrases[3]
    r4 = get_stressed_syllable(p4)
    results.append(r4)
    print(f"Phrase: '{p4}'")
    print(f"Rule: Mobile stress verb with no particles (Group 3), so stress is on the final syllable.")
    print(f"Result: The stressed syllable is number {r4}\n")

    # Phrase 5
    p5 = phrases[4]
    r5 = get_stressed_syllable(p5)
    results.append(r5)
    print(f"Phrase: '{p5}'")
    print(f"Rule: Verb has prefix 'vy-' (Group 1), so stress is on 'vy'.")
    print(f"Result: The stressed syllable is number {r5}\n")

    # Phrase 6
    p6 = phrases[5]
    r6 = get_stressed_syllable(p6)
    results.append(r6)
    print(f"Phrase: '{p6}'")
    print(f"Rule: Mobile stress verb with particle 'i' (Group 3), so stress is on the final syllable.")
    print(f"Result: The stressed syllable is number {r6}\n")

    print("-" * 50)
    print("The individual numbers for the final sequence are:")
    for num in results:
        print(num, end=" ")
    
    final_answer = "".join(map(str, results))
    print(f"\n\nThe final sequence of digits is: {final_answer}")
    
solve()