import itertools

# A small list of common English words for demonstration purposes
common_words = {
    "passes", "to", "the", "reception", "chamber"
}

# Scrambled sentence
scrambled_sentence = "peasss ot the ipcetnoer cebahrm"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)