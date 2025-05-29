# A small list of common English words for demonstration purposes
common_words = [
    "obtaining", "a", "copy", "upon", "request", "of", "the", "work", "in", "its", "original", "plain"
]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check each word in the common words list
    for word in common_words:
        if sorted(word) == sorted_scrambled:
            # Return the word with the original capitalization
            return word.capitalize() if scrambled_word[0].isupper() else word
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "ogtiinanb a cypo opun reqeust fo teh wokr ni sti lriginao Plani"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)