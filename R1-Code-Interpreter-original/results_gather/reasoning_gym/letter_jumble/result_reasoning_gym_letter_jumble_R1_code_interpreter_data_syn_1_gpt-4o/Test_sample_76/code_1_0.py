import itertools

# List of scrambled words
scrambled_words = ["si", "ilstl", "ni", "the", "ufrtue", "tbu", "ni"]

# A small set of common English words
valid_words = {"is", "still", "in", "the", "future", "but"}

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in valid_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)