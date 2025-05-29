# List of valid English words for this problem
valid_words = {"GAVE", "DAVE", "DOVE", "DOSE"}

# Sequence of transformations
word_sequence = ["GAVE", "DAVE", "DOVE", "DOSE"]

# Check if all words in the sequence are valid
all_valid = all(word in valid_words for word in word_sequence)

print(all_valid)