import itertools

# List of scrambled words
scrambled_words = ["yna", "omney", "pdia", "rof", "a", "wkro", "ro", "a", "retleaemcnp", "pocy", "fi"]

# Predefined list of common words
common_words = {
    "any", "money", "paid", "for", "a", "work", "or", "a", "replacement", "copy", "if"
}

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)