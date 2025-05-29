import itertools

# A small list of common English words for demonstration purposes
common_words = {
    "If", "you", "are", "outside", "the", "United", "States", "check", "laws", "of", "your", "country", "in", "addition", "to", "terms"
}

# List of scrambled words
scrambled_words = ["fI", "uoy", "era", "odesiut", "eht", "Utined", "Statse", "ccehk", "hte", "alws", "fo", "rouy", "tcunory", "ni", "andidiot", "ot", "hte", "tsrme", "fo"]

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
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)