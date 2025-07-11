def find_unstressed_words():
    """
    This script finds all words in a given Russian sentence that are not
    accented on the first syllable and have more than one syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # Russian vowels, both lowercase and uppercase. 'ё' is always stressed.
    vowels = "аеёиоуыэюяАЕЁИОУЫЭЮЯ"

    # A dictionary mapping lowercase words to the 0-based index of their stressed vowel.
    # This information is typically found in accentological dictionaries.
    stress_data = {
        # word: index_of_stressed_vowel
        "шашлык": 5,   # шашлЫк
        "запах": 1,    # зАпах
        "горелым": 3,  # горЕлым
        "вскоре": 2,   # вскОре
        "прибежал": 6, # прибежАл
        "сосед": 3,   # сосЕд
        "видно": 1,   # вИдно
        "дошёл": 3,    # дошЁл
        "его": 2,      # егО
        "дома": 1      # дОма (genitive case from "дом")
    }

    # Clean the text by replacing punctuation with spaces, then split into words.
    # This approach preserves the original capitalization of words.
    cleaned_text = text.replace(',', ' ').replace('-', ' ').replace('.', ' ')
    words = cleaned_text.split()

    result_words = []
    for word in words:
        if not word:
            continue
        
        lower_word = word.lower()

        # Find the indices of all vowels in the lowercase version of the word.
        vowel_indices = [i for i, char in enumerate(lower_word) if char in vowels]

        # 1. Exclude one-syllable words.
        if len(vowel_indices) <= 1:
            continue

        # 2. Exclude words where the stress is on the first syllable.
        if lower_word in stress_data:
            stressed_vowel_index = stress_data[lower_word]
            
            # The first syllable corresponds to the first vowel in the word.
            first_vowel_index = vowel_indices[0]

            if stressed_vowel_index != first_vowel_index:
                # The word is multi-syllable and not stressed on the first syllable.
                result_words.append(word)

    print(", ".join(result_words))

find_unstressed_words()