import itertools
import nltk
from nltk.corpus import words

# Ensure the words corpus is downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["dear", "rchAdbali", "keep", "no", "nSice", "yrstyedae", "atnhks", "ot", "oyu", "rheet", "si", "a", "gain", "fo", "sucrbbisesr", "rM"]

# Convert to lowercase for matching
scrambled_words_lower = [word.lower() for word in scrambled_words]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to find unscrambled word
def unscramble(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in valid_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words_lower]

# Reconstruct the sentence with the original capitalization and punctuation
unscrambled_sentence = " ".join(unscrambled_words)
print(unscrambled_sentence)