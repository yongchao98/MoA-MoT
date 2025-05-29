import itertools

# List of scrambled words
scrambled_words = [
    "rotfh", "ni", "tish", "egerenatm", "uoy", "mtus", "nboati", 
    "msrpsoeiin", "ni", "wigtnir", "ormf", "teh", "ePjorct", 
    "nureetGgb", "Leritary", "hrAcevi", "aoundonFti"
]

# A small dictionary of relevant English words
english_words = {
    "forth", "in", "this", "agreement", "you", "must", "obtain", 
    "permission", "in", "writing", "from", "the", "Project", 
    "Gutenberg", "Literary", "Archive", "Foundation"
}

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in english_words:
            return candidate
    return word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)