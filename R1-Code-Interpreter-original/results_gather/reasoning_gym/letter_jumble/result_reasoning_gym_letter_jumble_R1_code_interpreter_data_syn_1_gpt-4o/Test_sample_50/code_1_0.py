# A small list of common English words for demonstration
common_words = [
    "is", "was", "The", "Foundation", "business", "office"
]

# Scrambled sentence
scrambled_sentence = "s wlas The Foodnitaun s nbsisuse ofefci"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check against each word in the common words list
    for word in common_words:
        if sorted_scrambled == sorted(word.lower()):
            # Preserve the original capitalization
            if scrambled_word.istitle():
                return word.title()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)