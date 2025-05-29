import itertools
import nltk
from nltk.corpus import words

# Ensure the words corpus is downloaded
nltk.download('words')

# List of English words
word_list = set(words.words())

# Scrambled sentence
scrambled_sentence = "hte clal fo itearpemvi yudt gliht heatedr eh gpnissr fmro ihs ebd adn ertens"

# Function to unscramble a single word
def unscramble_word(scrambled_word):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate.lower() in word_list:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Output the result
print(unscrambled_sentence)