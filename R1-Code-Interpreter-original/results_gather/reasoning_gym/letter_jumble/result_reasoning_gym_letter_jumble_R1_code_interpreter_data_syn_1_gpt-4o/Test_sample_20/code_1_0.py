# Small dictionary of possible words
dictionary = [
    "regulating", "charities", "and", "charitable", "donations", "in", "all", "states", "of", "the", "United"
]

# List of scrambled words
scrambled_words = ["regulitang", "iharitces", "dna", "charitable", "donatnois", "ni", "all", "ttates", "fo", "eht", "Utined"]

# Function to unscramble a word
def unscramble(scrambled_word, dictionary):
    scrambled_sorted = sorted(scrambled_word)
    for word in dictionary:
        if sorted(word) == scrambled_sorted:
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word, dictionary) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)