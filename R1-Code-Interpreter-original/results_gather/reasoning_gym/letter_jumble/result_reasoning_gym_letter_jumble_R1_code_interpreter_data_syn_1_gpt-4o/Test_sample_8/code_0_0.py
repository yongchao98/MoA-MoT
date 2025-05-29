import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["iMehalc", "S", "arHt", "wsa", "eht", "iigortaron", "fo", "eht", "jrcPeto", "Gutnrbege", "ecntopc", "fo", "a", "iralbry", "fo", "eneriloctc", "skrow"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(word):
    # Generate all possible permutations of the word
    permutations = set(''.join(p) for p in itertools.permutations(word.lower()))
    # Find the correct word in the dictionary
    for perm in permutations:
        if perm in english_words:
            # Return the word with the original capitalization
            return perm.capitalize() if word[0].isupper() else perm
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)