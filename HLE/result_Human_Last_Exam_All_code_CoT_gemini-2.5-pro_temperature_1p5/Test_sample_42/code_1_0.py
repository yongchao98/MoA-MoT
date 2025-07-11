def solve():
    """
    Solves the word pattern puzzle by defining and applying a two-part rule.
    """
    vowels = "aeiou"
    forbidden_letters = ['m', 'p']

    def follows_pattern(word):
        """
        Checks if a word follows the identified pattern.
        Pattern:
        1. Must be a 3-letter word with a Consonant-Vowel-Consonant (CVC) structure.
        2. Must not contain any of the forbidden letters ('m', 'p').
        """
        # Rule 1 Check: CVC structure
        if len(word) != 3:
            return False, "Word is not 3 letters long."
        if word[0] in vowels or word[2] in vowels or word[1] not in vowels:
            return False, "Word does not follow the Consonant-Vowel-Consonant structure."

        # Rule 2 Check: Forbidden letters
        for letter in word:
            if letter in forbidden_letters:
                return False, f"Word contains a forbidden letter ('{letter}')."

        return True, "Word follows the pattern."

    answer_choices = {
        "A": "leg",
        "B": "dam",
        "C": "rat",
        "D": "car",
        "E": "bin",
    }

    print("Analyzing the answer choices...")
    non_conforming_word = None
    
    for key, word in answer_choices.items():
        does_follow, reason = follows_pattern(word)
        if does_follow:
            print(f"{key}. '{word}': Follows the pattern.")
        else:
            print(f"{key}. '{word}': Does NOT follow the pattern. Reason: {reason}")
            non_conforming_word = key

    print(f"\nThe word that does not follow the pattern is '{answer_choices[non_conforming_word]}'.")


solve()