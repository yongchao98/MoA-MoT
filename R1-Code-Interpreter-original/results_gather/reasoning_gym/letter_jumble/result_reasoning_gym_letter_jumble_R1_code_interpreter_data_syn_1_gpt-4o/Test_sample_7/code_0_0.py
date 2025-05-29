import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["sa", "amn", "himself", "rtecleicity", "aws", "knwon", "reays", "aog", "nad"]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate.lower() in valid_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)