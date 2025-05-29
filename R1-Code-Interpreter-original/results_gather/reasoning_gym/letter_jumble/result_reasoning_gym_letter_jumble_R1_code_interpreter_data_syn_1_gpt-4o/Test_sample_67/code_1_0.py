# A small list of common words for demonstration purposes
common_words = [
    "pointed", "to", "the", "sky", "which", "was", "of", "a", "pure", "blue", 
    "Yes", "muttered", "Mr", "Smith", "a"
]

# List of scrambled words
scrambled_words = ["dpinteo", "ot", "eht", "sk", "whihc", "wsa", "fo", "a", "puer", "bule", "eYs", "ttmreued", "rM", "himtS", "a"]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check each word in the common words list
    for word in common_words:
        if sorted(word.lower()) == sorted_scrambled:
            # Return the word with the original capitalization
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

print(sentence)