import re

def find_words_with_non_first_syllable_stress():
    """
    Finds and prints Russian words from a given text that have more than one
    syllable and where the stress is not on the first syllable.
    """
    # The input text in Russian
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A dictionary to map lowercase words to their stressed versions.
    # This is the most reliable method for determining stress in a fixed text.
    # The combining acute accent '́' follows the stressed vowel.
    # The letter 'ё' is inherently stressed.
    stress_map = {
        'шашлык': 'шашлы́к',
        'запах': 'за́пах',
        'горелым': 'горе́лым',
        'вскоре': 'вско́ре',
        'прибежал': 'прибежа́л',
        'сосед': 'сосе́д',
        'видно': 'ви́дно',
        'дошёл': 'дошёл',
        'его': 'его́',
        'дома': 'до́ма'
    }

    vowels = "аеёиоуыэюя"
    
    # Extract all sequences of Cyrillic letters, in order of appearance
    words = re.findall(r'[а-яА-ЯёЁ]+', text.lower())

    result_list = []
    for word in words:
        # Count syllables by counting vowels.
        syllable_count = sum(1 for char in word if char in vowels)

        # Rule: Exclude one-syllable words.
        if syllable_count <= 1:
            continue

        # Find the first vowel's position.
        first_vowel_pos = -1
        for i, char in enumerate(word):
            if char in vowels:
                first_vowel_pos = i
                break
        
        # Determine the stressed vowel's position.
        stressed_vowel_pos = -1
        if word in stress_map:
            stressed_word = stress_map[word]
            if '́' in stressed_word:
                # The stressed vowel is the character before the accent mark.
                stressed_vowel_pos = stressed_word.find('́') - 1
            elif 'ё' in stressed_word:
                # The letter 'ё' is always stressed.
                stressed_vowel_pos = stressed_word.find('ё')
        
        # If the stressed vowel's position is not the same as the first
        # vowel's position, the stress is not on the first syllable.
        if stressed_vowel_pos != -1 and stressed_vowel_pos != first_vowel_pos:
            result_list.append(word)

    # Print the final list, comma-separated.
    print(','.join(result_list))

find_words_with_non_first_syllable_stress()