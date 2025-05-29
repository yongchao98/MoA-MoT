import itertools

# List of scrambled words
scrambled_words = ["tiglkna", "ngonouripcn", "a", "anem", "hsi", "neam", "tizFr", "hTe", "lgidefthlu"]

# A small list of common English words for demonstration purposes
common_words = {
    "talking", "pronouncing", "a", "name", "his", "fritz", "the", "delightful"
}

# Function to unscramble a word
def unscramble(scrambled):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it is a valid English word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)