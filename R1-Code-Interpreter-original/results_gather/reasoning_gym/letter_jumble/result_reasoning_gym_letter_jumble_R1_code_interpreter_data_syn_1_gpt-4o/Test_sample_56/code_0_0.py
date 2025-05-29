import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["het", "hewol", "wrlod", "doucl", "eb", "eye", "sisnwsete", "fo", "het", "necfeorrmpa", "a", "petorrer", "mwehaelin", "eikl"]

# Get the list of valid English words
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in valid_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)