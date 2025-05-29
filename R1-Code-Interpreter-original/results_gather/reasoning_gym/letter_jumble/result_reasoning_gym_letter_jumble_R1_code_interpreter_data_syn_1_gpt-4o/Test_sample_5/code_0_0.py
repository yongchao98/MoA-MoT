import itertools

# A small dictionary of words for demonstration purposes
dictionary = {"stand", "in", "isolation", "where", "we"}

def unscramble_word(scrambled, dictionary):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        word = ''.join(perm)
        if word in dictionary:
            return word
    return scrambled  # Return the original if no match is found

# The scrambled sentence
scrambled_sentence = "nstaodino ni olcatosni eehrw ew"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)