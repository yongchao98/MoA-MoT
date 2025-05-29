# A small list of common English words for demonstration purposes
common_words = [
    "from", "several", "printed", "editions", "all", "of", "which", "are", 
    "confirmed", "as", "not", "protected", "by", "copyright", "in", "the", "U", "S"
]

# List of scrambled words
scrambled_words = [
    "mrof", "selerav", "drintep", "ediniots", "all", "fo", "wihch", "aer", 
    "confermid", "sa", "nto", "proteetcd", "yb", "copyrtghi", "ni", "hte", "U", "S"
]

# Function to unscramble a word
def unscramble_word(scrambled, word_list):
    scrambled_sorted = sorted(scrambled)
    for word in word_list:
        if sorted(word) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word, common_words) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)