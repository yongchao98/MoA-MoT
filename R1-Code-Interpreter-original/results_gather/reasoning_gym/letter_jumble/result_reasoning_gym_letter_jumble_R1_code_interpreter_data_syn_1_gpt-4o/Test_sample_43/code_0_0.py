import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "rotfh", "ni", "tish", "egerenatm", "uoy", "mtus", "nboati", 
    "msrpsoeiin", "ni", "wigtnir", "ormf", "teh", "ePjorct", 
    "nureetGgb", "Leritary", "hrAcevi", "aoundonFti"
]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)