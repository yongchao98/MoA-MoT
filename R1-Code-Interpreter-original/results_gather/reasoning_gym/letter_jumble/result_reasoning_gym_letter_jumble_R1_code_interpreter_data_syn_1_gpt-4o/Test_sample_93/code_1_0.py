import itertools

# A small dictionary of common words for demonstration purposes
common_words = {
    "clock", "answered", "Mr", "Smith", "and", "with", "that", "the", "day", "is"
}

# Scrambled sentence
scrambled_sentence = "ckocl arswdnee rM Stimh adn iwth ahtt hte ady s"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word):
        possible_word = ''.join(perm)
        if possible_word in common_words:
            return possible_word
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)