# A small list of common English words for demonstration purposes
common_words = {
    "and", "physically", "forces", "depend", "on", "the", "mode", "of", 
    "vibration", "of", "the", "either", "particles", "which", "is", 
    "for", "each", "specifically", "different"
}

# Scrambled sentence
scrambled_sentence = "dna psihacyl efrcos penedd no eht omde fo viirnbtoa fo hte erheitc ptrlicsea whcih si fro aech sceciiafplly difrnfeet"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble_word(scrambled):
    for word in common_words:
        if sorted(word) == sorted(scrambled):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)