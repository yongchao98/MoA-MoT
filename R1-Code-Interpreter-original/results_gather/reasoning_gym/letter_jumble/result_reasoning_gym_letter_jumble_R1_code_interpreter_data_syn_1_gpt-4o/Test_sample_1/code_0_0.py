import itertools

# A simple dictionary of words for demonstration purposes
dictionary = ["Why", "not", "Mr", "Smith", "advanced", "to", "the", "young"]

# Scrambled words
scrambled_words = ["hWy", "nto", "rM", "imSth", "adncvade", "ot", "hte", "yuong"]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled.lower()) == sorted(word.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

print(unscrambled_sentence)