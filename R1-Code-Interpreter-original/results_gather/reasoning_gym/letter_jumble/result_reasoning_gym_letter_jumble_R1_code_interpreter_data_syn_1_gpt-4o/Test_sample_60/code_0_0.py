import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["eht", "uorf", "prcjetos", "fyboralav", "ceevired", "eht", "first", "aws", "ttha", "fo", "a", "gonuy", "amn", "shwoe", "badro", "feaehodr"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in english_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)