import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["tkhin", "thne", "reeakdrm", "rM", "hSimt", "lsuyngim", "that", "fi", "one", "docul"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        # Check if the candidate is a valid English word
        if candidate.lower() in english_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)