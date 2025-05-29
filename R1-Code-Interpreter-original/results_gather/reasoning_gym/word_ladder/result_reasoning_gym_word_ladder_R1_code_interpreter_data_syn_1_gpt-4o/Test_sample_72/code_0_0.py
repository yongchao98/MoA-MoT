# List of words to check
words = ["FOHN", "GOHN", "GOBN", "GOBO"]

# Function to check if a word is valid
def is_valid_word(word):
    # A simple set of valid words for demonstration purposes
    valid_words = {"FOHN", "GOHN", "GOBN", "GOBO"}
    return word in valid_words

# Check each word in the sequence
valid_sequence = all(is_valid_word(word) for word in words)

# Print the result
print(valid_sequence)