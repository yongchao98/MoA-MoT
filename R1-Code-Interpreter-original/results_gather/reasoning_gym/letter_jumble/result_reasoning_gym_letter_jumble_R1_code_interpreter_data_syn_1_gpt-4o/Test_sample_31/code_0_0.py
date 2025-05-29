import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "ym", "omclrucatau", "wrsko", "ta", "aNagria", "Fro", "rM", "htimS", 
    "ahs", "boiatend", "a", "lease", "fo", "hte", "arget", "sflla", "fo"
]

# Convert the list of words to a set for faster lookup
word_set = set(words.words())

# Function to unscramble a word
def unscramble(scrambled):
    scrambled_lower = scrambled.lower()
    for perm in itertools.permutations(scrambled_lower):
        candidate = ''.join(perm)
        if candidate in word_set:
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)