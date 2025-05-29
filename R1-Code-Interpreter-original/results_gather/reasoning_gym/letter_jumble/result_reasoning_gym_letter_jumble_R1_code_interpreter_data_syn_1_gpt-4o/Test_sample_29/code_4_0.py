# List of scrambled words
scrambled_words = ["taht", "I", "uowld", "od", "inythang", "os", "pnratuiotic", "eW", "ta", "eomh", "tnihk", "taht", "yuor", "pelpoe", "muts", "won", "eb", "setad"]

# A basic list of common English words
common_words = [
    "that", "I", "would", "do", "anything", "so", "patriotic", "We", "at", "home", "think", "your", "people", "must", "now", "be", "stead"
]

# Function to unscramble a word
def unscramble(word, common_words):
    for common_word in common_words:
        if sorted(word.lower()) == sorted(common_word.lower()):
            return common_word
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word, common_words) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)