import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["oFrm", "hte", "ady", "ehnw", "shete", "wto", "inntmusetrs", "were", "cvrtoedin", "si", "ot", "eb", "dated", "eht"]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word.lower()))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in valid_words:
            # Return the word with the original capitalization
            if scrambled_word[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    # If no valid word is found, return the original scrambled word
    return scrambled_word

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)