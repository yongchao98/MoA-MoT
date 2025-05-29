import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["eDfects", "suhc", "sa", "btu", "nto", "ldmitei", "ot", "itcomplene", "iuaccnrate", "ro"]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word.lower()):
        candidate = ''.join(perm)
        if candidate in valid_words:
            # Preserve the original capitalization
            if scrambled_word[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)