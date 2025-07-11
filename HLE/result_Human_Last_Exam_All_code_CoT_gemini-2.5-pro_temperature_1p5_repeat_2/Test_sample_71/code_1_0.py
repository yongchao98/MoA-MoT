import re

def solve_old_russian_stress():
    """
    Solves the Old Russian stress puzzle by applying a deduced set of phonological rules.
    The code explains the logic and calculates the stressed syllable for each target phrase.
    """
    
    print("This script determines the stressed syllable in Old Russian phrases based on a set of derived rules.")
    print("The rules are applied in the following order of priority:")
    print("1. 'vy-' Prefix Rule: The prefix 'vy-' is always stressed.")
    print("2. Negation Rule: If 'ne' is present with endings '-lo' or '-li', 'ne' is stressed.")
    print("3. Verb Class Rule: Otherwise, stress depends on the verb's class (Root-stressed vs. Ending-stressed).")
    print("-" * 20)

    target_phrases = {
        "i ne znali": "Rule 2 applies ('ne' + '-li'). Stress on 'ne'.",
        "i povelo že": "Rule 3, Class II ('ved-'). Stress on the last syllable 'že'.",
        "ne vymyla že": "Rule 1 applies. Stress on prefix 'vy-'.",
        "ponesla": "Rule 3, Class II ('nes-'). Stress on the ending '-la'.",
        "vyvela že": "Rule 1 applies. Stress on prefix 'vy-'.",
        "i unesli": "Rule 3, Class II ('nes-'). Stress on the ending '-li'."
    }
    
    results = []
    
    for i, (phrase, explanation) in enumerate(target_phrases.items()):
        vowels = re.findall(r'[aeiouy]', phrase)
        syllable_count = len(vowels)
        
        # Determine stressed syllable based on the pre-analyzed rules
        if i == 0: # i ne znali
            stress_pos = phrase.find('ne')
            syllables_before = len(re.findall(r'[aeiouy]', phrase[:stress_pos]))
            result = syllables_before + 1
        elif i == 1: # i povelo že
            result = syllable_count
        elif i == 2: # ne vymyla že
            stress_pos = phrase.find('vy')
            syllables_before = len(re.findall(r'[aeiouy]', phrase[:stress_pos]))
            result = syllables_before + 1
        elif i == 3: # ponesla
            result = syllable_count
        elif i == 4: # vyvela že
            stress_pos = phrase.find('vy')
            syllables_before = len(re.findall(r'[aeiouy]', phrase[:stress_pos]))
            result = syllables_before + 1
        elif i == 5: # i unesli
            result = syllable_count

        results.append(result)
        
        print(f"Phrase: '{phrase}'")
        print(f"Analysis: {explanation}")
        print(f"Result: The stressed syllable is {result}")
        print("-" * 20)

    final_answer = "".join(map(str, results))
    print("Final sequence of stressed syllable numbers:")
    # Per instruction "output each number in the final equation", this loop prints each digit.
    for num in final_answer:
        print(num, end="")
    print() # for a final newline

solve_old_russian_stress()
<<<252314>>>