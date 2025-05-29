import itertools

# List of scrambled words
scrambled_words = ["etxenivse", "ni", "eth", "wdrol", "cegishtrnt"]

# A small list of common English words for demonstration
common_words = {"extensive", "in", "the", "world", "stretching"}

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)