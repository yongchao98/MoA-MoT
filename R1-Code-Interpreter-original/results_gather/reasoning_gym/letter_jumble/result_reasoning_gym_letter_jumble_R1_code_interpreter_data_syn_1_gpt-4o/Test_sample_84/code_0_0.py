import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "arYe", "asw", "tirsf", "eshlibpdu", "ni", "eth", "rebuaryF", "p", 
    "tI", "asw", "epsibluhd", "ni", "Fnerca", "het", "tnxe", "year", 
    "thghouAl", "pdlbsihue", "nrued", "eth"
]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled):
    scrambled_lower = scrambled.lower()
    for word in english_words:
        if sorted(scrambled_lower) == sorted(word):
            return word
    return scrambled

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
sentence = " ".join(unscrambled_words)

# Print the unscrambled sentence
print(sentence)