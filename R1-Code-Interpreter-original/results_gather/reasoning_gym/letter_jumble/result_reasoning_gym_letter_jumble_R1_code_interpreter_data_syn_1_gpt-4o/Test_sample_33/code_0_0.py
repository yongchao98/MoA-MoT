import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["fsee", "ot", "mete", "nad", "eekp", "pu", "wtih", "shete", "risueremqnte", "eW", "od", "nto", "sclioti", "donotiasn", "ni", "iocaltons"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)