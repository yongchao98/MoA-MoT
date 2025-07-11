def solve_stress():
    """
    This script determines the stressed syllable in six Old Russian phrases
    based on a derived set of rules.
    """

    # --- Rule Set ---
    # 1. Verbs are classified as Group A (stem-stressed) or Group B (end-stressed).
    #    - Group A: 'zna', 'myl', and any verb with prefix 'vy-'.
    #    - Group B: 'nes', 'vel'.
    # 2. Stress placement rules based on particles:
    #    - No particles, or with 'i' only: Stress is on the verb's inherent position (stem for A, end for B).
    #    - With 'ne' only: Stress shifts to 'ne' for Group B verbs.
    #    - With 'že' only: Stress shifts to 'že' for Group B verbs.
    #    - Prefixed 'vy-' verbs are always stressed on 'vy'.
    #    - Particle 'i' seems to cancel the stress-shifting effect of 'ne' and 'že'.

    vowels = "aeiouy"

    def find_vowel_position(phrase, stressed_vowel_char_index):
        """Counts vowels from the left to find the position of the stressed one."""
        vowel_count = 0
        for i, char in enumerate(phrase):
            if char in vowels:
                vowel_count += 1
            if i == stressed_vowel_char_index:
                return vowel_count
        return -1

    results = []
    print("Determining the stressed syllable for each phrase:")
    print("-" * 50)

    # 1. i ne znali
    phrase = "i ne znali"
    # Root 'zna' is Group A. Particle 'i' cancels 'ne' effect. Rule: stem stress.
    # Stressed vowel is 'a' in 'znali'.
    stressed_char_index = phrase.find('a')
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the root is 'zna' (Group A). With 'i ne', stress remains on the stem 'zna'.")
    print(f"The stressed vowel is the 3rd vowel ('a'). Stressed syllable: {stress_pos}\n")

    # 2. i povelo že
    phrase = "i povelo že"
    # Root 'vel' is Group B. Particle 'i' cancels 'že' effect. Rule: end stress.
    # Stressed vowel is the ending 'o' in 'povel-o'.
    stressed_char_index = phrase.rfind('o', 0, phrase.find('že'))
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the root is 'vel' (Group B). With 'i...že', stress falls on the ending 'o'.")
    print(f"The stressed vowel is the 4th vowel ('o'). Stressed syllable: {stress_pos}\n")

    # 3. ne vymyla že
    phrase = "ne vymyla že"
    # Prefix 'vy-' is always stressed (Group A behavior).
    # Stressed vowel is 'y' in 'vy-'.
    stressed_char_index = phrase.find('y')
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the prefix 'vy-' is always stressed.")
    print(f"The stressed vowel is the 2nd vowel ('y'). Stressed syllable: {stress_pos}\n")

    # 4. ponesla
    phrase = "ponesla"
    # Root 'nes' is Group B. No particles. Rule: end stress.
    # Stressed vowel is the ending 'a' in 'ponesl-a'.
    stressed_char_index = phrase.find('a')
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the root is 'nes' (Group B). With no particles, stress is on the ending 'a'.")
    print(f"The stressed vowel is the 3rd vowel ('a'). Stressed syllable: {stress_pos}\n")
    
    # 5. vyvela že
    phrase = "vyvela že"
    # Prefix 'vy-' is always stressed (Group A behavior).
    # Stressed vowel is 'y' in 'vy-'.
    stressed_char_index = phrase.find('y')
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the prefix 'vy-' is always stressed.")
    print(f"The stressed vowel is the 1st vowel ('y'). Stressed syllable: {stress_pos}\n")

    # 6. i unesli
    phrase = "i unesli"
    # Root 'nes' is Group B. Particle 'i'. Rule: end stress.
    # Stressed vowel is the ending 'i' in 'unesl-i'.
    stressed_char_index = phrase.rfind('i')
    stress_pos = find_vowel_position(phrase, stressed_char_index)
    results.append(stress_pos)
    print(f"In '{phrase}', the root is 'nes' (Group B). With 'i', stress is on the ending 'i'.")
    print(f"The stressed vowel is the 4th vowel ('i'). Stressed syllable: {stress_pos}\n")

    print("-" * 50)
    final_answer = "".join(map(str, results))
    print(f"The final combined number is: {final_answer}")

solve_stress()