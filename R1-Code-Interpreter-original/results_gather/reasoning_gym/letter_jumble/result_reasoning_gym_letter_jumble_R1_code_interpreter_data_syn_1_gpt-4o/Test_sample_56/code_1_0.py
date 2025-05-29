import itertools

# List of scrambled words
scrambled_words = ["het", "hewol", "wrlod", "doucl", "eb", "eye", "sisnwsete", "fo", "het", "necfeorrmpa", "a", "petorrer", "mwehaelin", "eikl"]

# A small set of common English words for demonstration
common_words = {
    "the", "whole", "world", "could", "be", "eye", "witness", "of", "the", "performance", "a", "reporter", "meanwhile", "like"
}

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)