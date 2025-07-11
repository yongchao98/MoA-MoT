def solve():
    """
    This script identifies the word that does not follow a specific pattern.
    The pattern is:
    1. The word must have a Consonant-Vowel-Consonant (CVC) structure.
    2. The word must not contain the letters 'm' or 'p'.
    The script checks each answer choice against this pattern.
    """
    vowels = "aeiou"
    answer_choices = {
        "A": "leg",
        "B": "dam",
        "C": "rat",
        "D": "car",
        "E": "bin"
    }

    offending_word = None
    offending_key = None

    print("Analyzing answer choices based on the pattern (CVC structure, no 'm' or 'p'):\n")

    for key, word in answer_choices.items():
        # Check CVC structure
        is_cvc = (
            len(word) == 3 and
            word[0] not in vowels and
            word[1] in vowels and
            word[2] not in vowels
        )

        # Check for forbidden letters
        has_forbidden_letter = 'm' in word or 'p' in word

        # A word follows the pattern if it's CVC AND has no forbidden letters
        if is_cvc and not has_forbidden_letter:
            print(f"'{word}' ({key}): Follows the pattern.")
        else:
            offending_word = word
            offending_key = key
            reason = []
            if not is_cvc:
                reason.append("it does not have a Consonant-Vowel-Consonant structure")
            if has_forbidden_letter:
                reason.append(f"it contains a forbidden letter ('m' or 'p')")
            print(f"'{word}' ({key}): Does NOT follow the pattern because {' and '.join(reason)}.")

    if offending_word:
        print(f"\nThe word that does not follow the pattern is '{offending_word}'.")

solve()