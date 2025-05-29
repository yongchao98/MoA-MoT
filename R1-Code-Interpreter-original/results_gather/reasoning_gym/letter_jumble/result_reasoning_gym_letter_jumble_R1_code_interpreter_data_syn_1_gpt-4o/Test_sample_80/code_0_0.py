import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["tse", "orf", "shi", "reoirtncuesr", "ti", "swa", "poposrde", "ot", "rM", "hSmti", "ttah", "eh", "ohsudl", "rpmiet", "eth", "snecdo", "tpar", "fo"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)