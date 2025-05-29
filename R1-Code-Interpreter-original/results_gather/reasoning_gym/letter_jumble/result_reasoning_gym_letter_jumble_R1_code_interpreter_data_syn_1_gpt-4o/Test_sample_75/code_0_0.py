import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "any", "part", "of", "this", "collection", "work", "without", "employment", "displaying", "the"
]

# Scrambled words
scrambled_words = ["ayn", "tpra", "fo", "iths", "ecritlonec", "work", "ihtwout", "ilyeonpntrm", "dgsiiaynlp", "teh"]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        # Find a word in the dictionary that matches the scrambled word
        for dict_word in dictionary:
            if sorted(word) == sorted(dict_word):
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
unscrambled_sentence_str = ' '.join(unscrambled_sentence)

print(unscrambled_sentence_str)