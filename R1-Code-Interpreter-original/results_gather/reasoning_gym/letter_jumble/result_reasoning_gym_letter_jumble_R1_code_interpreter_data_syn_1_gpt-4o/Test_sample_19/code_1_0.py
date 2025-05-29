import itertools

# List of scrambled words
scrambled_words = ["tkhin", "thne", "reeakdrm", "rM", "hSimt", "lsuyngim", "that", "fi", "one", "docul"]

# A small set of common English words for demonstration
common_words = {
    "think", "then", "remarked", "Mr", "Smith", "smiling", "that", "if", "one", "could"
}

# Function to unscramble a word
def unscramble(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        # Check if the candidate is a valid English word
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)