# List of scrambled words
scrambled_words = ["ayn", "tpra", "fo", "iths", "ecritlonec", "work", "ihtwout", "ilyeonpntrm", "dgsiiaynlp", "teh"]

# A small dictionary of common English words
common_words = {
    "any", "part", "of", "this", "collection", "work", "without", "employment", "displaying", "the"
}

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = ''.join(sorted(scrambled_word))
    # Check each word in the common words list
    for word in common_words:
        if ''.join(sorted(word)) == sorted_scrambled:
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)