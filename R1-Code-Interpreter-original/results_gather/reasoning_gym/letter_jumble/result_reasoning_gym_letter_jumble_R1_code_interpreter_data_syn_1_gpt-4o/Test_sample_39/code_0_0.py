import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# The scrambled sentence
scrambled_sentence = "Verne ordRatce s neot eht Year saw frits upblidhes ni teh rbFyuare p"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled_word):
    # Check all permutations of the scrambled word
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        # Check if the candidate is a valid English word
        if candidate.lower() in english_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)