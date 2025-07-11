import re
try:
    from russtress import Accent
except ImportError:
    print("Error: The 'russtress' library is not installed.")
    print("Please install it by running: pip install russtress")
    exit()

# The original text in Russian
text = "Шашлык запах горелым, и вскоре прибежал сосед - видно, запах дошёл и до его дома."

# A string containing all Russian vowels for syllable counting
vowels = "аеёиоуыэюя"

# Initialize the accent detector from the russtress library
accent = Accent()

# Clean the text: replace punctuation with spaces, then split into a list of words
# This preserves the original order.
cleaned_text = re.sub(r'[,\-.]', ' ', text)
words = cleaned_text.split()

# This list will store the words that meet our criteria
result_words = []

# Iterate through each word from the cleaned text
for word in words:
    # Use the lowercase version for analysis to ensure consistency
    word_lower = word.lower()
    
    # Count the number of syllables by counting the vowels in the word
    syllable_count = sum(1 for char in word_lower if char in vowels)
    
    # Condition 1: The word must have more than one syllable.
    if syllable_count > 1:
        
        # Find the index of the first vowel in the word. This marks the start of the first syllable.
        first_vowel_index = -1
        for i, char in enumerate(word_lower):
            if char in vowels:
                first_vowel_index = i
                break
        
        # Assume the stress is not on the first syllable until proven otherwise
        is_first_syllable_stressed = False
        
        # Use the library to place a stress mark (') in the word
        stressed_word = accent.put_stress(word_lower, stress_symbol="'")
        
        # Find the position of the stress mark
        stress_mark_index = stressed_word.find("'")
        
        # If a stress mark was placed by the library...
        if stress_mark_index != -1:
            # Check if the stress mark immediately follows the first vowel.
            # If so, the stress is on the first syllable.
            if stress_mark_index == first_vowel_index + 1:
                is_first_syllable_stressed = True
        else:
            # As a fallback, handle words with 'ё'. The letter 'ё' is always stressed.
            # If 'ё' is the first vowel of the word, then the stress is on the first syllable.
            if 'ё' in word_lower and word_lower.find('ё') == first_vowel_index:
                is_first_syllable_stressed = True
        
        # Condition 2: If the stress is NOT on the first syllable, add the original word to our result list.
        if not is_first_syllable_stressed:
            result_words.append(word)

# Print the final list of words, separated by commas, as requested.
print(",".join(result_words))