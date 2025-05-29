import itertools

# A simple dictionary of words for demonstration purposes
dictionary = {"in", "journalism", "well", "cash"}

def unscramble_word(scrambled_word, dictionary):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        word = ''.join(perm)
        if word.lower() in dictionary:
            return word
    return scrambled_word  # Return the original if no match is found

# The scrambled sentence
scrambled_sentence = "ni oiurnalmsj Well Cahs"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)