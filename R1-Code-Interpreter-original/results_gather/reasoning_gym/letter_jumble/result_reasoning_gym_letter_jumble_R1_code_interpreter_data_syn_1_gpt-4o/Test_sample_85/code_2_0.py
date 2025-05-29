import nltk
from nltk.corpus import words

# Download the word list if not already downloaded
nltk.download('words')

# Use the nltk words corpus as the dictionary
dictionary = set(words.words())

# Scrambled words
scrambled_words = ["fo", "teh", "eeoplp", "fo", "areGt", "Birtain", "hwo", "tlli", "nwo", "ahd"]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for word in dictionary:
        if sorted(word.lower()) == sorted(scrambled.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

print(unscrambled_sentence)