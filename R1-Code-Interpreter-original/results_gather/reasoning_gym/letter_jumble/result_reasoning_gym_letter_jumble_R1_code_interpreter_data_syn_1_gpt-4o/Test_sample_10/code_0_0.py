import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of English words
word_list = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in word_list:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "ew hsall eb ebla ot puodrce"

# Unscramble each word in the sentence
unscrambled_words = [unscramble(word) for word in scrambled_sentence.split()]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)