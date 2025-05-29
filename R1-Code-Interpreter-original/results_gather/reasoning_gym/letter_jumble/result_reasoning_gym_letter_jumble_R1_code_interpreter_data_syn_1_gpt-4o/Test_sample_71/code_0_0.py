import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "tyeh", "are", "erev", "aseridgnsd", "ot", "het", "nearicAm", 
    "nmveonergt", "ivna", "atoilmcspn", "A", "imcpaagn", "stgiaan", 
    "eth", "anaaenxnito", "thatt", "has"
]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled):
    for perm in itertools.permutations(scrambled):
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)