import itertools

# List of scrambled words
scrambled_words = ["whti", "a", "miosnwe", "Wrehe", "haev"]

# A small set of common English words for demonstration
valid_words = {"with", "a", "someone", "where", "have"}

# Function to unscramble a word
def unscramble_word(scrambled):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in valid_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)