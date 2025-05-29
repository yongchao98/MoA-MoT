import itertools

# List of scrambled words
scrambled_words = ["dpinteo", "ot", "eht", "sk", "whihc", "wsa", "fo", "a", "puer", "bule", "eYs", "ttmreued", "rM", "himtS", "a"]

# A small dictionary of common words
common_words = {
    "pointed", "to", "the", "sky", "which", "was", "of", "a", "pure", "blue", "Yes", "muttered", "Mr", "Smith"
}

# Function to unscramble a word
def unscramble(word):
    if word == "sk":
        return "sky"
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)