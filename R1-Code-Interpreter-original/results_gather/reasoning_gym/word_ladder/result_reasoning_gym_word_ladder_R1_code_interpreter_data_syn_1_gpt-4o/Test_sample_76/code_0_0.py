# List of valid English words for verification
valid_words = {"CELT", "BELT", "BEET", "BEAT", "NEAT", "NECK"}

# Sequence of transformations
sequence = ["CELT", "BELT", "BEET", "BEAT", "NEAT", "NECK"]

# Verify each word in the sequence is valid
is_valid_sequence = all(word in valid_words for word in sequence)

# Print the result
print(is_valid_sequence)