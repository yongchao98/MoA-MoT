def solve_stress_puzzle():
    """
    Solves the Old Russian stress puzzle by applying a derived set of rules.
    """

    phrases = [
        "i ne znali",
        "i povelo że",
        "ne vymyla że",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    results = []

    # A helper function to count the number of vowels up to a specific character index
    def count_vowels_until(text, char_index):
        vowel_count = 0
        for i, char in enumerate(text):
            if char in "aeiouy":
                vowel_count += 1
            if i == char_index:
                return vowel_count
        return vowel_count

    # Rule application for each phrase

    # 1. i ne znali: Fixed root stress on 'znal'. Stress on 'a'.
    phrase = phrases[0]
    stress_pos = phrase.find('zna') + 2 # Position of 'a' in 'zna'
    results.append(str(count_vowels_until(phrase, stress_pos)))

    # 2. i povelo że: Mobile stress, '-lo' form with 'że'. Stress on 'że'.
    phrase = phrases[1]
    vowels = [c for c in phrase if c in "aeiouy"]
    results.append(str(len(vowels))) # 'że' is the last syllable

    # 3. ne vymyla że: Prefix 'vy-' is always stressed. Stress on 'y'.
    phrase = phrases[2]
    stress_pos = phrase.find('vy') + 1 # Position of 'y' in 'vy'
    results.append(str(count_vowels_until(phrase, stress_pos)))
    
    # 4. ponesla: Mobile stress, '-la' form. Default stress on the ending 'la'. Stress on 'a'.
    phrase = phrases[3]
    stress_pos = phrase.find('la') + 1 # Position of 'a' in 'la'
    results.append(str(count_vowels_until(phrase, stress_pos)))

    # 5. vyvela że: Prefix 'vy-' is always stressed. Stress on 'y'.
    phrase = phrases[4]
    stress_pos = phrase.find('vy') + 1 # Position of 'y' in 'vy'
    results.append(str(count_vowels_until(phrase, stress_pos)))

    # 6. i unesli: Mobile stress, '-li' form, no retracting particles. Default stress on ending 'li'. Stress on 'i'.
    phrase = phrases[5]
    # Find the last 'i' in the phrase which corresponds to the ending.
    stress_pos = phrase.rfind('i')
    results.append(str(count_vowels_until(phrase, stress_pos)))
    
    print("".join(results))

solve_stress_puzzle()
<<<352314>>>