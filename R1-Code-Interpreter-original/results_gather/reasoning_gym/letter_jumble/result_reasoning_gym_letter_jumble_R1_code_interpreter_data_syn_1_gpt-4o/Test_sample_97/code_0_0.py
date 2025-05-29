import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = ["fI", "uoy", "era", "odesiut", "eht", "Utined", "Statse", "ccehk", "hte", "alws", "fo", "rouy", "tcunory", "ni", "andidiot", "ot", "hte", "tsrme", "fo"]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid English word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in english_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)