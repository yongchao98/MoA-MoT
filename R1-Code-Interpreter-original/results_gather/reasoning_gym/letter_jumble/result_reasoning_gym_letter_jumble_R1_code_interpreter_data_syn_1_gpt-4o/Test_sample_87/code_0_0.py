import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of English words
word_list = set(words.words())

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all possible permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate.lower() in word_list:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "ogtiinanb a cypo opun reqeust fo teh wokr ni sti lriginao Plani"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)