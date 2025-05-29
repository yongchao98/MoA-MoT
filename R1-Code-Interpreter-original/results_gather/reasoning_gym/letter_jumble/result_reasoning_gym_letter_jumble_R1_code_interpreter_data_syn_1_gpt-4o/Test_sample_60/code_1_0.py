# A small list of common English words for demonstration purposes
common_words = {
    "the", "four", "projects", "favorably", "received", "first", "was", "that", "of", "a", "young", "man", "whose", "board", "forehead"
}

# List of scrambled words
scrambled_words = ["eht", "uorf", "prcjetos", "fyboralav", "ceevired", "eht", "first", "aws", "ttha", "fo", "a", "gonuy", "amn", "shwoe", "badro", "feaehodr"]

# Function to unscramble a word
def unscramble(scrambled_word):
    from itertools import permutations
    for perm in permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)