# List of words to verify
words = ["TOWS", "SOWS", "SOIS", "SOIL"]

# A simple set of valid words for demonstration purposes
valid_words = {"TOWS", "SOWS", "SOIL"}

# Check if each word in the sequence is valid
valid_sequence = all(word in valid_words for word in words)

# Print the result
print(valid_sequence)