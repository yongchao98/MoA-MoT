import itertools

# A simple dictionary of words for demonstration purposes
dictionary = {"The", "first", "thing", "that"}

# Function to unscramble a word
def unscramble(word):
    # Generate all permutations of the word
    permutations = set(''.join(p) for p in itertools.permutations(word))
    # Find the correct word in the dictionary
    for perm in permutations:
        if perm in dictionary:
            return perm
    return word  # Return the original word if no match is found

# Scrambled sentence
scrambled_sentence = "heT isrft nhigt hatt"

# Unscramble each word in the sentence
unscrambled_words = [unscramble(word) for word in scrambled_sentence.split()]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)