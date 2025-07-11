import re

def find_words_with_non_initial_stress():
    """
    From the given Russian text, this function lists all words (excluding one-syllable words),
    comma-separated and in the order of appearance, that do not have the accent
    on the first syllable.
    """
    # The input text in Russian.
    text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."
    
    # A dictionary containing pre-analyzed stress information for multi-syllable words.
    # The key is the lowercase word, and the value is the 0-indexed position of the stressed vowel.
    # For example, in 'шашлык' (vowels 'а', 'ы'), the stress is on 'ы' (шашлЫк),
    # which is the second vowel, so its index is 1.
    stress_data = {
        "шашлык":   1,  # шаш-лЫк
        "запах":    0,  # зА-пах
        "горелым":  1,  # го-рЕ-лым
        "вскоре":   0,  # вскО-ре
        "прибежал": 2,  # при-бе-жАл
        "сосед":   1,  # со-сЕд
        "видно":   0,  # вИд-но
        "дошёл":    1,  # до-шЁл
        "его":      1,  # е-гО
        "дома":     0,  # дО-ма
    }
    
    # Define Russian vowels to count syllables.
    vowels = "аеёиоуыэюя"
    
    # Use regex to find all sequences of Cyrillic letters, preserving order.
    words = re.findall(r'[а-яА-ЯёЁ]+', text)
    
    result_words = []
    # Iterate through the words found in the text.
    for word in words:
        # Count vowels to approximate syllable count.
        syllable_count = sum(1 for char in word.lower() if char in vowels)
        
        # Exclude one-syllable words.
        if syllable_count <= 1:
            continue
            
        word_lower = word.lower()
        
        # Check if the word is in our data and its stress is not on the first syllable (index > 0).
        if word_lower in stress_data and stress_data[word_lower] > 0:
            result_words.append(word)
            
    # Print the final list, with words separated by a comma and a space.
    print(", ".join(result_words))

find_words_with_non_initial_stress()
<<<Шашлык, горелым, прибежал, сосед, дошёл, его>>>