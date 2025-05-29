import itertools

# List of scrambled words
scrambled_words = ["F", "Seom", "atsest", "od", "not", "awlol", "iidseacmlrs", "fo", "nearitc"]

# A small set of common English words for demonstration
common_words = {
    "a", "some", "states", "do", "not", "allow", "disclaimers", "of", "certain", "f"
}

# Function to unscramble a word
def unscramble(word):
    if word.lower() in common_words:
        return word
    for perm in itertools.permutations(word.lower()):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)