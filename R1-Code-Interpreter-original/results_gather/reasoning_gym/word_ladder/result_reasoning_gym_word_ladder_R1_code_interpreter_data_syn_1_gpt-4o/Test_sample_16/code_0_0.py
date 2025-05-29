# List of words to verify
words = ["DICE", "TICE", "TINE", "TEEN"]

# A simple set of valid words for demonstration purposes
valid_words = {"DICE", "TICE", "TINE", "TEEN"}

# Check if all words in the sequence are valid
all_valid = all(word in valid_words for word in words)

print(all_valid)