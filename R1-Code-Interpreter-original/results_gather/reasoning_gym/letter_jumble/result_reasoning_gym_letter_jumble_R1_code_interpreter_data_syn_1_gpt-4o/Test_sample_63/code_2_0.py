import itertools

# List of scrambled words
scrambled_words = [
    "os", "atht", "ti", "yam", "ees", "untesrdand", "adn", "mbrereem", "teh",
    "reasnos", "which", "eteirmend", "hte", "pelionrysat", "ttah", "atsc",
    "tJsu", "tyuds", "uyolsfre", "sa"
]

# A basic list of common English words with corrections
common_words = {
    "so", "that", "it", "may", "see", "understand", "and", "remember", "the",
    "reasons", "which", "determined", "the", "personality", "that", "acts",
    "Just", "study", "yourself", "as"
}

# Function to unscramble a single word
def unscramble_word(scrambled):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)