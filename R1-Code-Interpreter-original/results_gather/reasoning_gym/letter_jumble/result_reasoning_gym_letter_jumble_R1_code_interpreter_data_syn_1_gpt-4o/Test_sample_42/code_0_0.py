import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "yb", "iosnelptCro", "emti", "ouy", "neam", "sYe", "oodG", "yb", "tnhe", 
    "rfo", "a", "tetlli", "hielw", "idas", "rM", "mtSih", "sa", "eh", "ervdees"
]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for perm in itertools.permutations(scrambled):
        candidate = ''.join(perm)
        if candidate.lower() in dictionary:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word, english_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

print(sentence)