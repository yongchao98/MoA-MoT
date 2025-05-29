# List of words to check
words = ["LEAP", "SEAP", "SEEP", "SEED"]

# Function to check if a word is valid
def is_valid_word(word):
    # A simple set of valid words for demonstration purposes
    valid_words = {"LEAP", "SEAP", "SEEP", "SEED"}
    return word in valid_words

# Check each word in the sequence
valid_sequence = all(is_valid_word(word) for word in words)

# Print the result
print(valid_sequence)