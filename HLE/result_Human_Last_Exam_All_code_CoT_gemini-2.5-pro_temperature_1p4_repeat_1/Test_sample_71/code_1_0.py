def solve_stress():
    """
    This script determines the stressed syllable in Old Russian phrases based on a set of derived rules.
    """

    phrases_to_solve = {
        "i ne znali": "Rule B: i + ne -> stress on last syllable of the verb.",
        "i povelo že": "Rule D: že suffix, prefixed verb -> stress on verb root 've'.",
        "ne vymyla že": "Rule C: ne + že -> stress on first syllable of the verb 'vy'.",
        "ponesla": "Rule E: prefixed verb -> stress on verb root 'nes'.",
        "vyvela že": "Rule D: že suffix, prefixed verb -> stress on verb prefix 'vy'.",
        "i unesli": "Rule E: i + prefixed verb -> stress on verb root 'nes'."
    }
    
    vowels = "aeiouy"
    final_answer_digits = []

    # 1. i ne znali
    # Rule B: Stress on the last syllable of the verb. The last vowel is the 'i' in 'znali'.
    phrase = "i ne znali"
    syllable_count = 0
    for char in phrase:
        if char in vowels:
            syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")

    # 2. i povelo že
    # Rule D: Stress on the verb root 've'.
    phrase = "i povelo že"
    target_syllable_part = "ve"
    syllable_count = 0
    # Find the start of the target syllable and count vowels up to that point.
    pos = phrase.find(target_syllable_part)
    # Count vowels in the substring before the target part
    for char in phrase[:pos]:
        if char in vowels:
            syllable_count += 1
    # Add one for the target syllable itself
    syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")
    
    # 3. ne vymyla že
    # Rule C: Stress on the first syllable of the verb, 'vy'.
    phrase = "ne vymyla že"
    target_syllable_part = "vy"
    syllable_count = 0
    pos = phrase.find(target_syllable_part)
    for char in phrase[:pos]:
        if char in vowels:
            syllable_count += 1
    syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")
    
    # 4. ponesla
    # Rule E: Stress on the verb root, 'nes'.
    phrase = "ponesla"
    target_syllable_part = "nes"
    syllable_count = 0
    pos = phrase.find(target_syllable_part)
    for char in phrase[:pos]:
        if char in vowels:
            syllable_count += 1
    syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")

    # 5. vyvela že
    # Rule D: Stress on the verb prefix, 'vy'.
    phrase = "vyvela že"
    target_syllable_part = "vy"
    syllable_count = 0
    pos = phrase.find(target_syllable_part)
    for char in phrase[:pos]:
        if char in vowels:
            syllable_count += 1
    syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")

    # 6. i unesli
    # Rule E: Stress on the verb root, 'nes'.
    phrase = "i unesli"
    target_syllable_part = "nes"
    syllable_count = 0
    pos = phrase.find(target_syllable_part)
    for char in phrase[:pos]:
        if char in vowels:
            syllable_count += 1
    syllable_count += 1
    results = {"phrase": phrase, "pos": syllable_count}
    final_answer_digits.append(str(results['pos']))
    print(f"'{results['phrase']}': The stressed syllable is {results['pos']}")

    final_answer = "".join(final_answer_digits)
    print(f"\nThe final answer is: {final_answer}")
    
solve_stress()