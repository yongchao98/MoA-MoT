import re

def find_words_with_non_initial_stress():
    """
    Analyzes a Russian sentence to find all multi-syllable words
    where the stress is not on the first syllable.
    """
    # The original Russian text
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # Russian vowels are used to count syllables
    vowels = "аоуыэяёюие"

    # A map of words from the text to their stressed versions.
    # The character 'ё' is always stressed, so words with it are handled separately.
    # The stress is marked by the combining acute accent character (U+0301).
    stress_map = {
        'шашлык': 'шашлы́к',
        'запах': 'за́пах',
        'горелым': 'горе́лым',
        'вскоре': 'вско́ре',
        'прибежал': 'прибежа́л',
        'сосед': 'сосе́д',
        'видно': 'ви́дно',
        'его': 'его́',
        'дома': 'дома́',  # In "до его дома" (genitive case), stress is on the 'а'.
    }

    # Tokenize the text: split by any sequence of non-Russian letters.
    # This handles spaces, commas, hyphens, and periods correctly.
    # We also filter out empty strings that can result from the split.
    words = [word for word in re.split(r'[^а-яё]+', text.lower()) if word]

    # This list will hold the words that match the criteria.
    result_words = []

    # Process each word in its order of appearance.
    for word in words:
        # Count syllables by counting the vowels.
        syllable_count = sum(1 for char in word if char in vowels)

        # Exclude one-syllable words as per the requirement.
        if syllable_count <= 1:
            continue

        # Find the character index of the first vowel in the word.
        first_vowel_index = -1
        for i, char in enumerate(word):
            if char in vowels:
                first_vowel_index = i
                break
        
        # Determine the index of the stressed vowel.
        stressed_vowel_index = -1
        if 'ё' in word:
            # If 'ё' is present, it's the stressed vowel.
            stressed_vowel_index = word.find('ё')
        elif word in stress_map:
            # Otherwise, use the map to find the stressed vowel.
            stressed_word = stress_map[word]
            # The stressed vowel's index in the original word is one less than
            # the index of the combining accent mark.
            stressed_vowel_index = stressed_word.find('\u0301') - 1

        # Check if the stressed vowel is the first vowel of the word.
        # If it's not the first, the stress is not on the first syllable.
        if stressed_vowel_index != first_vowel_index:
            result_words.append(word)

    # Print the final list of words, comma-separated.
    print(", ".join(result_words))

find_words_with_non_initial_stress()