import itertools

# List of scrambled words
scrambled_words = [
    "remvoing", "na", "eneirt", "ciyt", "lal", "ta", "oecn", "fmor", "eno", 
    "plaec", "ot", "naother", "Hsi", "speaicl", "orpject", "hda", "ot", "od", "iwth", "teh"
]

# A small set of common English words for demonstration purposes
common_words = {
    "removing", "an", "entire", "city", "all", "at", "once", "from", "one", 
    "place", "to", "another", "his", "special", "project", "had", "to", "do", "with", "the"
}

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            # Preserve the original capitalization
            if word.istitle():
                return candidate.capitalize()
            else:
                return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)