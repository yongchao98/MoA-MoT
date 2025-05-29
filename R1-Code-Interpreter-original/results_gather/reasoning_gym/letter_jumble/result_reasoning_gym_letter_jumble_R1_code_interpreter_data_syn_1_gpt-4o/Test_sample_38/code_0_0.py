import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of English words
word_list = set(words.words())

# Scrambled sentence
scrambled_sentence = "ntcoain Dtfeecs scuh sa ubt ton"

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        possible_word = ''.join(perm)
        if possible_word.lower() in word_list:
            return possible_word
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)