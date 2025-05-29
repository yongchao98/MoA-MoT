from itertools import permutations

# A simple dictionary of words for demonstration purposes
dictionary = {"definitely", "settled", "That", "great"}

def unscramble_word(scrambled_word):
    # Generate all possible permutations of the scrambled word
    possible_words = {''.join(p) for p in permutations(scrambled_word)}
    # Find the intersection with the dictionary
    matches = possible_words.intersection(dictionary)
    return matches.pop() if matches else scrambled_word

# Scrambled sentence
scrambled_sentence = ["definytiveli", "settdel", "Ttah", "greta"]

# Unscramble each word
unscrambled_sentence = [unscramble_word(word) for word in scrambled_sentence]

# Join the unscrambled words into a sentence
result_sentence = ' '.join(unscrambled_sentence)

print(result_sentence)