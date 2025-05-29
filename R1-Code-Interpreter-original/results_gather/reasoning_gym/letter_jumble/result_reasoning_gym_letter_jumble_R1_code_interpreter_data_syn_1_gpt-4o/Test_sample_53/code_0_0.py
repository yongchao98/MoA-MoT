import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["tiglkna", "ngonouripcn", "a", "anem", "hsi", "neam", "tizFr", "hTe", "lgidefthlu"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it is a valid English word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)