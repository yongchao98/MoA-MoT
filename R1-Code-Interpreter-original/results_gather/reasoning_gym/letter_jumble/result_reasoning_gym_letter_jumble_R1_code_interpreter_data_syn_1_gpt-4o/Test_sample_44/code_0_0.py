import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_sentence = "tinUed Sattes nad oyu rae dotecal ni teh dniUet States ew od ton cmali"
scrambled_words = scrambled_sentence.split()

# Get the list of valid words from nltk
valid_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled):
        candidate = ''.join(perm)
        # Check if the candidate is a valid word
        if candidate.lower() in valid_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)