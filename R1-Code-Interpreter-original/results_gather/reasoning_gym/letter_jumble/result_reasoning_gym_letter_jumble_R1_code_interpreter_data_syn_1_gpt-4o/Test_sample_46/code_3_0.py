# List of common English words
common_words = [
    "strong", "to", "recalls", "One", "readily", "understands", "how", "a", "man"
]

# List of scrambled words
scrambled_words = ["grastn", "ot", "calrels", "enO", "realidy", "undersnatds", "how", "a", "nam"]

# Function to unscramble a word
def unscramble(scrambled, word_list):
    scrambled_sorted = sorted(scrambled.lower())
    for word in word_list:
        if sorted(word.lower()) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word, common_words) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)