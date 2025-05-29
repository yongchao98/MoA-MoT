import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "os", "atht", "ti", "yam", "ees", "untesrdand", "adn", "mbrereem", "teh", 
    "reasnos", "which", "eteirmend", "hte", "pelionrysat", "ttah", "atsc", 
    "tJsu", "tyuds", "uyolsfre", "sa"
]

# Get the list of valid words from nltk corpus
valid_words = set(words.words())

# Function to unscramble a single word
def unscramble_word(scrambled):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in valid_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)