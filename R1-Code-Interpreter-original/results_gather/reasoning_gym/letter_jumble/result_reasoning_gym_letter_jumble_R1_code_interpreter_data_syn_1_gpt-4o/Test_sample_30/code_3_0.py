import itertools

# A small list of common English words for demonstration purposes
common_words = {
    "this", "for", "like", "all", "wealthy", "folk", "in", "our"
}

# List of scrambled words
scrambled_words = ["dhsi", "Fro", "lkie", "lla", "aewlthy", "oflk", "ni", "uor"]

# Function to unscramble a word
def unscramble(scrambled):
    scrambled_lower = scrambled.lower()
    for perm in itertools.permutations(scrambled_lower):
        candidate = ''.join(perm)
        if candidate in common_words:
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)