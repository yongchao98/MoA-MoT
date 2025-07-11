def solve_old_russian_stress():
    """
    This script determines the stressed syllable in Old Russian phrases based on a derived set of rules.
    """
    # Define morpheme properties based on the analysis of examples
    root_types = {'zna': 'A', 'my': 'A', 'nes': 'B', 've': 'B'}
    prefixes = ['vy', 'po', 'u']
    roots = ['zna', 'my', 'nes', 've']
    endings = ['la', 'lo', 'li']

    # The phrases to analyze
    phrases_to_analyze = [
        "i ne znali",
        "i povelo že",
        "ne vymyla že",
        "ponesla",
        "vyvela že",
        "i unesli"
    ]

    final_answer_digits = []

    print("Determining the stressed syllable for each phrase:")

    for phrase_str in phrases_to_analyze:
        # 1. Parse the phrase into its component morphemes (syllables)
        components = []
        temp_phrase = phrase_str

        # Handle particles at the start
        if temp_phrase.startswith("i "):
            components.append("i")
            temp_phrase = temp_phrase[2:]
        if temp_phrase.startswith("ne "):
            components.append("ne")
            temp_phrase = temp_phrase[2:]

        # Handle particle at the end
        has_ze_particle = False
        if temp_phrase.endswith(" že"):
            has_ze_particle = True
            temp_phrase = temp_phrase[:-3]

        # Handle prefix
        found_prefix = None
        for p in prefixes:
            if temp_phrase.startswith(p):
                found_prefix = p
                components.append(p)
                temp_phrase = temp_phrase[len(p):]
                break
        
        # Handle root and ending
        found_root = None
        found_ending = None
        # The remainder of the phrase is root + ending
        for r in roots:
            if temp_phrase.startswith(r):
                found_root = r
                components.append(r)
                ending_part = temp_phrase[len(r):]
                if ending_part in endings:
                    found_ending = ending_part
                    components.append(found_ending)
                break
        
        if has_ze_particle:
            components.append("že")

        # 2. Apply the hierarchical stress rules
        stress_syllable_index = -1

        # Rule 1: Prefix 'vy-' is always stressed
        if 'vy' in components:
            stress_syllable_index = components.index('vy')
        # Rule 2: 'ne' is stressed unless preceded by 'i'
        elif 'ne' in components and (components.index('ne') == 0 or components[components.index('ne')-1] != 'i'):
            stress_syllable_index = components.index('ne')
        # Rule 3: Stress depends on root type (default case)
        else:
            if root_types.get(found_root) == 'A':  # Type A: root-stressed
                stress_syllable_index = components.index(found_root)
            elif root_types.get(found_root) == 'B':  # Type B: ending-stressed
                stress_syllable_index = components.index(found_ending)
        
        # Result is the 1-based syllable number
        result = stress_syllable_index + 1
        final_answer_digits.append(str(result))
        
        # Output the result for the current phrase
        print(f"'{phrase_str}': {result}")

    # Combine the digits for the final answer
    final_answer = "".join(final_answer_digits)
    print(f"\nThe final combined answer is: {final_answer}")

# Execute the function
solve_old_russian_stress()