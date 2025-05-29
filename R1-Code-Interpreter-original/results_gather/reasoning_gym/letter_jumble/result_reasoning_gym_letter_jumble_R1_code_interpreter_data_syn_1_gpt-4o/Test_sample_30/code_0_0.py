import itertools
import nltk
from nltk.corpus import words

# Ensure the words corpus is downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["dhsi", "Fro", "lkie", "lla", "aewlthy", "oflk", "ni", "uor"]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled):
    scrambled_lower = scrambled.lower()
    for perm in itertools.permutations(scrambled_lower):
        candidate = ''.join(perm)
        if candidate in valid_words:
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)