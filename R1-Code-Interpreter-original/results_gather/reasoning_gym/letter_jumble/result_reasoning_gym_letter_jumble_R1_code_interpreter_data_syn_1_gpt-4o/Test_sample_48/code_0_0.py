import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["wlil", "eb", "eknild", "ot", "eht", "rPojtce", "Guerngetb", "Lsceine", "fro", "lal", "osrkw"]

# Get the list of valid words from the nltk corpus
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        # Check if the candidate is a valid word
        if candidate.lower() in valid_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)