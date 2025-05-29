# List of words to verify
words = ["TOWS", "SOWS", "SOLS", "SOIL"]

# A more comprehensive set of valid words
valid_words = {"TOWS", "SOWS", "SOLS", "SOIL"}

# Check if each word in the sequence is valid
valid_sequence = all(word in valid_words for word in words)

# Print the result
print(valid_sequence)