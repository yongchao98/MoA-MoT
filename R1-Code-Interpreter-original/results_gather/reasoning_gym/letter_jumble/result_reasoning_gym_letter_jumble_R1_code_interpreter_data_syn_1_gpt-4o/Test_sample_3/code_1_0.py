import itertools

# A small sample list of common English words for demonstration
common_words = {
    "Mr", "Smith", "continues", "his", "round", "and", "enters", "the", 
    "reporters", "hall", "Here", "reporters", "in", "their", "respective", 
    "places", "facing", "an", "equal", "number"
}

# List of scrambled words
scrambled_words = [
    "rM", "Smtih", "inntsouec", "ihs", "dounr", "nad", "tneers", "eht", 
    "rpotreers", "ahll", "reHe", "rreoptesr", "ni", "theri", "eesprcteiv", 
    "slpcea", "gacnif", "na", "eauql", "buenmr"
]

# Function to unscramble a word
def unscramble(scrambled):
    for perm in itertools.permutations(scrambled):
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
sentence = " ".join(unscrambled_words)

print(sentence)