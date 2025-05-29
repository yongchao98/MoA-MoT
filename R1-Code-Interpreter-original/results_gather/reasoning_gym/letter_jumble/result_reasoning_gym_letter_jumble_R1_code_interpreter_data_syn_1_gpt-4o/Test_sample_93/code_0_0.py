import itertools

# A small dictionary of words for demonstration purposes
dictionary = set([
    "clock", "answered", "Mr", "Smith", "and", "with", "that", "the", "day", "is"
])

# The scrambled sentence
scrambled_sentence = "ckocl arswdnee rM Stimh adn iwth ahtt hte ady s"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Generate all possible permutations of the scrambled word
    permutations = itertools.permutations(scrambled_word)
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        possible_word = ''.join(perm)
        if possible_word in dictionary:
            return possible_word
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)