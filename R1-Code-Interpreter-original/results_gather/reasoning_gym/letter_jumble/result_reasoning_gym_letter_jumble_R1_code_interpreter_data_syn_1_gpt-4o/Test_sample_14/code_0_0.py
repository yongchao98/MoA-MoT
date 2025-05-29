import itertools
import nltk
from nltk.corpus import words

# Ensure the words corpus is downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ['ni', 'eht', 'nognimr', 'a', 'atelb']

# Get the list of valid English words
english_words = set(words.words())

# Function to unscramble a single word
def unscramble(word):
    for perm in itertools.permutations(word):
        candidate = ''.join(perm)
        if candidate in english_words:
            return candidate
    return word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)