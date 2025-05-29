# Predefined list of common English words
common_words = [
    "Why", "not", "Mr", "Smith", "advanced", "to", "the", "young"
]

# List of scrambled words
scrambled_words = ["hWy", "nto", "rM", "imSth", "adncvade", "ot", "hte", "yuong"]

# Function to find unscrambled word
def unscramble(scrambled, word_list):
    scrambled_lower = scrambled.lower()
    for word in word_list:
        if sorted(scrambled_lower) == sorted(word.lower()):
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return word.capitalize()
            else:
                return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word, common_words) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)