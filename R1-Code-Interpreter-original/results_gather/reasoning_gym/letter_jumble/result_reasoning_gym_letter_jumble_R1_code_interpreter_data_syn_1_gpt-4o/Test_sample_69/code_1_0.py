# A small list of common English words for demonstration purposes
common_words = [
    "hypocrity", "aerostyle", "special", "rules", "set", "forth", "in", "the", "general"
]

# List of scrambled words
scrambled_words = ["htcporgiy", "rsyeloati", "aepclSi", "rslue", "tse", "hfrto", "ni", "eth", "Gaenerl"]

# Function to unscramble a word
def unscramble(scrambled_word, word_list):
    scrambled_word_lower = scrambled_word.lower()
    for word in word_list:
        if sorted(scrambled_word_lower) == sorted(word):
            # Preserve the original capitalization
            if scrambled_word[0].isupper():
                return word.capitalize()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word, common_words) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)