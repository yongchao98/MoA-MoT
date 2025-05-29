# List of words in the sequence
words = ["REEK", "PEEK", "PUKE", "PUCK"]

# A simple set of valid English words for demonstration purposes
valid_words = {"REEK", "PEEK", "PUKE", "PUCK"}

# Check if all words in the sequence are valid
all_valid = all(word in valid_words for word in words)

# Print the result
print(all_valid)