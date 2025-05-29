import itertools

# A small list of common English words for demonstration purposes
common_words = {
    "department", "he", "found", "the", "operators", "testing", "with", "folded", "arms", "at", "their"
}

# List of scrambled words
scrambled_words = ["redapmtent", "eh", "ofudn", "eht", "oprtoaers", "istnigt", "hwti", "ofledd", "samr", "ta", "thire"]

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)