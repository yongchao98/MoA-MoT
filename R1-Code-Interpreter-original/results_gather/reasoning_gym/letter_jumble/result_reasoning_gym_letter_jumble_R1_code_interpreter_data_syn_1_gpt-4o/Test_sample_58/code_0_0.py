import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["F", "Seom", "atsest", "od", "not", "awlol", "iidseacmlrs", "fo", "nearitc"]

# Get the list of valid English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    if word.lower() in english_words:
        return word
    for perm in itertools.permutations(word.lower()):
        candidate = ''.join(perm)
        if candidate in english_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)