import itertools

# List of scrambled words
scrambled_words = ["taht", "I", "uowld", "od", "inythang", "os", "pnratuiotic", "eW", "ta", "eomh", "tnihk", "taht", "yuor", "pelpoe", "muts", "won", "eb", "setad"]

# A small list of common English words for demonstration purposes
common_words = {
    "that", "I", "would", "do", "anything", "so", "patriotic", "We", "at", "home", "think", "your", "people", "must", "now", "be", "stead"
}

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)