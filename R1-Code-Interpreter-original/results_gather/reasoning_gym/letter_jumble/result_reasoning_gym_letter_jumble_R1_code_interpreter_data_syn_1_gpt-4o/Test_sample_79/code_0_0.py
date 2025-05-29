import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already available
nltk.download('words')

# List of scrambled words
scrambled_words = ["teeeohpnl", "I", "was", "wtih", "rD", "Wlkisni", "hA", "the", "dogo"]

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    if word.lower() in valid_words:
        return word
    for perm in itertools.permutations(word):
        perm_word = ''.join(perm)
        if perm_word.lower() in valid_words:
            return perm_word
    return word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)