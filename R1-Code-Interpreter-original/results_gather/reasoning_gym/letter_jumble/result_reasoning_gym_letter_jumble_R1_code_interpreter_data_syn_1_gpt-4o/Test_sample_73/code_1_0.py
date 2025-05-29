# A small sample list of common English words for demonstration
common_words = {
    "the", "call", "of", "imperative", "duty", "light", "hearted", "he", "springs", "from", "his", "bed", "and", "enters"
}

# Scrambled sentence
scrambled_sentence = "hte clal fo itearpemvi yudt gliht heatedr eh gpnissr fmro ihs ebd adn ertens"

# Function to unscramble a single word
def unscramble_word(scrambled_word):
    from itertools import permutations
    for perm in permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Output the result
print(unscrambled_sentence)