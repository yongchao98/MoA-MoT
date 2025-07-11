import re

def solve():
    """
    Finds and prints Russian words from a given text that have more than one syllable
    and are not stressed on the first syllable.
    """
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

    # A knowledge base mapping lowercase words to their stressed versions.
    # A capital letter indicates the stressed vowel. 'ё' is always stressed.
    stress_data = {
        "шашлык": "шашлЫк",
        "запах": "зАпах",
        "горелым": "горЕлым",
        "вскоре": "вскОре",
        "прибежал": "прибежАл",
        "сосед": "сосЕд",
        "видно": "вИдно",
        "дошёл": "дошЁл",
        "его": "егО",
        "дома": "дОма",
    }
    
    # Russian vowels for syllable counting.
    vowels = "аеёиоуыэюя"

    # Extract all words from the text in order of appearance.
    words = re.findall(r'\b[а-яА-ЯёЁ]+\b', text)

    result_words = []
    for word in words:
        word_lower = word.lower()

        # Rule 1: Exclude one-syllable words.
        syllable_count = sum(1 for char in word_lower if char in vowels)
        if syllable_count <= 1:
            continue

        # Rule 2: Exclude words with stress on the first syllable.
        if word_lower in stress_data:
            stressed_form = stress_data[word_lower]

            # Find the index of the first vowel.
            first_vowel_index = -1
            for i, char in enumerate(word_lower):
                if char in vowels:
                    first_vowel_index = i
                    break
            
            # Find the index of the stressed vowel.
            stressed_vowel_index = -1
            if 'ё' in stressed_form.lower():
                stressed_vowel_index = stressed_form.lower().find('ё')
            else:
                for i, char in enumerate(stressed_form):
                    if char.isupper():
                        stressed_vowel_index = i
                        break

            # If the stressed vowel is not the first vowel, add the original word to the result.
            if stressed_vowel_index != first_vowel_index:
                result_words.append(word)

    print(",".join(result_words))

solve()
<<<Шашлык,горелым,прибежал,сосед,дошёл,его>>>