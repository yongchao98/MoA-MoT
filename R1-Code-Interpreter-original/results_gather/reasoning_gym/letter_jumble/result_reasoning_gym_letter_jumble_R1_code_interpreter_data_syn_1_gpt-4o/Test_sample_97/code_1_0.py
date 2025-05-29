# A small list of common English words for demonstration purposes
common_words = {
    "if", "you", "are", "outside", "the", "united", "states", "check", "laws", "of", "your", "country", "in", "addition", "to", "terms"
}

# List of scrambled words
scrambled_words = ["fI", "uoy", "era", "odesiut", "eht", "Utined", "Statse", "ccehk", "hte", "alws", "fo", "rouy", "tcunory", "ni", "andidiot", "ot", "hte", "tsrme", "fo"]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = ''.join(sorted(scrambled_word.lower()))
    # Check against the sorted letters of each word in the common words list
    for word in common_words:
        if sorted_scrambled == ''.join(sorted(word)):
            # Preserve the original capitalization
            return word.capitalize() if scrambled_word[0].isupper() else word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)