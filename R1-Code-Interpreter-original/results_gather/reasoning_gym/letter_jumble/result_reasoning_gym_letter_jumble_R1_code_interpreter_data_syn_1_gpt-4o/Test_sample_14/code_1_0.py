import itertools

# List of scrambled words
scrambled_words = ['ni', 'eht', 'nognimr', 'a', 'atelb']

# A small set of common English words for demonstration
common_english_words = {'in', 'the', 'morning', 'a', 'table'}

# Function to unscramble a single word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in common_english_words:
            return candidate
    return word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)