import itertools

# A small list of common English words for demonstration purposes
common_words = {"the", "opposite", "side", "who", "knows"}

def unscramble_word(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

def unscramble_sentence(sentence):
    scrambled_words = sentence.split()
    unscrambled_words = [unscramble_word(word) for word in scrambled_words]
    return ' '.join(unscrambled_words)

# The scrambled sentence
scrambled_sentence = "teh osppeiot dsie hwo wnsko"
unscrambled_sentence = unscramble_sentence(scrambled_sentence)

print(unscrambled_sentence)