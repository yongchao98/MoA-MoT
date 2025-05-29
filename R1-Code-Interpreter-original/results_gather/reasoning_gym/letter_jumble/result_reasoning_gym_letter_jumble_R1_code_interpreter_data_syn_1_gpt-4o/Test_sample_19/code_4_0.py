# A small list of common English words for demonstration purposes
common_words = [
    "think", "then", "remarked", "Mr", "Smith", "smilingly", "that", "if", "one", "could"
]

# List of scrambled words
scrambled_words = ["tkhin", "thne", "reeakdrm", "rM", "hSimt", "lsuyngim", "that", "fi", "one", "docul"]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check each word in the common words list
    for word in common_words:
        if sorted(word.lower()) == sorted_scrambled:
            # Return the word with the original capitalization
            return word if scrambled_word[0].islower() else word.capitalize()
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)