# List of words to verify
words = ["LASS", "PASS", "PARS", "PART", "PORT"]

# Function to check if a word is valid
def is_valid_word(word):
    # A simple set of valid words for demonstration purposes
    valid_words = {"LASS", "PASS", "PARS", "PART", "PORT"}
    return word in valid_words

# Verify each word in the sequence
valid_sequence = all(is_valid_word(word) for word in words)

# Print the result
print(valid_sequence)