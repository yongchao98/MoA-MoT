import itertools

# List of scrambled words
scrambled_words = ["omon", "ouabt", "eacf", "hA", "teher", "s", "oimenshtg"]

# A small set of common words for demonstration purposes
common_words = {"moon", "about", "face", "a", "there", "is", "something", "ha"}

# Function to unscramble a word
def unscramble(scrambled):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in common_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)