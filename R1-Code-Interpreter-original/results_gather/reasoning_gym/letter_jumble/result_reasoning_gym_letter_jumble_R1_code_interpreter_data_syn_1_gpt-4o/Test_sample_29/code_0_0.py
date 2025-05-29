import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["taht", "I", "uowld", "od", "inythang", "os", "pnratuiotic", "eW", "ta", "eomh", "tnihk", "taht", "yuor", "pelpoe", "muts", "won", "eb", "setad"]

# Get the list of valid English words
valid_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in valid_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)