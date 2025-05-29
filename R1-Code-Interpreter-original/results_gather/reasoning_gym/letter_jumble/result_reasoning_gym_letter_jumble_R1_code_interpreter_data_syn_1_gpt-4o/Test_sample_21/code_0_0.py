import itertools

# A simple dictionary of words for demonstration purposes
dictionary = ["forth", "in", "this", "agreement", "violates", "the", "law", "of", "the"]

# Scrambled words
scrambled_words = ["tohfr", "ni", "thsi", "mganreeet", "volitaes", "teh", "lwa", "fo", "hte"]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        for dict_word in dictionary:
            if sorted(word) == sorted(dict_word):
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

print(final_sentence)