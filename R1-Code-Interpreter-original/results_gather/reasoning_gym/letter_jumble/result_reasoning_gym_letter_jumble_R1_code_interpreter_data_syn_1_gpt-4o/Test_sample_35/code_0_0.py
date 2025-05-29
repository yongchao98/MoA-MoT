import itertools

# A simple dictionary of words for demonstration purposes
dictionary = {"copy", "it", "give", "away"}

def unscramble_word(scrambled_word, dictionary):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check which permutation is a valid word in the dictionary
    for perm in permutations:
        word = ''.join(perm)
        if word in dictionary:
            return word
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "cypo ti giev ti away"
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)