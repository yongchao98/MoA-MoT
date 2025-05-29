import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "broke", "in", "Mr", "Smith", "turning", "to", "a", "second", "scientist", 
    "Optical", "science", "defective"
]

# Scrambled words
scrambled_words = [
    "borke", "ni", "rM", "mSith", "ntruing", "ot", "a", "sdcone", 
    "ncseitist", "tpOilac", "cniesce", "defcetive"
]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = sorted(word.lower())
        for dict_word in dictionary:
            if sorted_word == sorted(dict_word.lower()):
                # Preserve capitalization
                if word[0].isupper():
                    unscrambled_sentence.append(dict_word.capitalize())
                else:
                    unscrambled_sentence.append(dict_word)
                break
    return " ".join(unscrambled_sentence)

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)
print(unscrambled_sentence)