from itertools import permutations

# A small dictionary of words for demonstration purposes
dictionary = {"the", "opposite", "side", "who", "knows"}

def unscramble_word(scrambled, dictionary):
    # Generate all permutations of the scrambled word
    perms = {''.join(p) for p in permutations(scrambled)}
    # Find the intersection of permutations and the dictionary
    matches = perms.intersection(dictionary)
    return matches.pop() if matches else scrambled

# Scrambled sentence
scrambled_sentence = "teh osppeiot dsie hwo wnsko"
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)
print(unscrambled_sentence)