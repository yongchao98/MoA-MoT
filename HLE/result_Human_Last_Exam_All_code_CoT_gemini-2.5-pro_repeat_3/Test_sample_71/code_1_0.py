import re

def solve_stress():
    """
    This function determines the stressed syllable in Old Russian phrases based on a set of derived rules.
    """
    
    # Helper function to count syllables by finding vowels.
    # Returns the 1-based syllable number for a given character index.
    def get_syllable_number(phrase, char_index):
        vowels = "aeiouy"
        count = 0
        for i, char in enumerate(phrase):
            if char in vowels:
                count += 1
                if i == char_index:
                    return count
        return -1 # Should not happen if char_index is a vowel index

    # The phrases to be analyzed
    phrases = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]
    
    # Store the results
    results = []

    print("Determining stress for each phrase:")

    # 1. 'i ne znali'
    # Analysis: Root is '-zna-', which has fixed stress. Rule 2 applies.
    # The stress should be on the vowel 'a' of the root 'zna'.
    phrase = phrases[0]
    stress_index = phrase.find('zna') + 1 # Index of 'a' in 'zna'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"1. In 'i ne znali', the root '-zna-' has fixed stress. The stress falls on the root's vowel 'a', which is syllable number: {result}")

    # 2. 'i povelo že'
    # Analysis: Root is '-ve-', which has mobile stress. It's a neuter form with 'že'. Rule 3b applies.
    # The stress should be on the particle 'že'.
    phrase = phrases[1]
    stress_index = phrase.find('že') # Index of 'e' in 'že'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"2. In 'i povelo že', the root '-ve-' has mobile stress and 'že' is present. The stress falls on 'že', which is syllable number: {result}")

    # 3. 'ne vymyla že'
    # Analysis: The verb contains the prefix 'vy-'. Rule 1 applies.
    # The stress should be on the prefix 'vy-'.
    phrase = phrases[2]
    stress_index = phrase.find('vy') + 1 # Index of 'y' in 'vy'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"3. In 'ne vymyla že', the prefix 'vy-' is present. The stress falls on 'vy-', which is syllable number: {result}")

    # 4. 'ponesla'
    # Analysis: Root is '-nes-', which has mobile stress. It's a feminine form ending in '-a'. Rule 3a applies.
    # The stress should be on the final ending '-a'.
    phrase = phrases[3]
    stress_index = phrase.rfind('a') # Index of the final 'a'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"4. In 'ponesla', the root '-nes-' has mobile stress and the form is feminine. The stress falls on the ending '-a', which is syllable number: {result}")

    # 5. 'vyvela že'
    # Analysis: The verb contains the prefix 'vy-'. Rule 1 applies.
    # The stress should be on the prefix 'vy-'.
    phrase = phrases[4]
    stress_index = phrase.find('vy') + 1 # Index of 'y' in 'vy'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"5. In 'vyvela že', the prefix 'vy-' is present. The stress falls on 'vy-', which is syllable number: {result}")

    # 6. 'i unesli'
    # Analysis: Root is '-nes-', which has mobile stress. It's a plural form without 'ne' or 'že'. Rule 3c applies.
    # The stress should be on the final ending '-i'.
    phrase = phrases[5]
    stress_index = phrase.rfind('i') # Index of the final 'i'
    result = get_syllable_number(phrase, stress_index)
    results.append(result)
    print(f"6. In 'i unesli', the root '-nes-' has mobile stress and no particles 'ne'/'že'. The stress falls on the ending '-li', which is syllable number: {result}")

    # Combine the results into a single string
    final_answer = "".join(map(str, results))
    print("\nThe final sequence of stressed syllable numbers is:")
    print(final_answer)

solve_stress()