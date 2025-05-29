import itertools
from nltk.corpus import words

# Load the list of English words
word_list = set(words.words())

# Scrambled sentence
scrambled_sentence = "dna psihacyl efrcos penedd no eht omde fo viirnbtoa fo hte erheitc ptrlicsea whcih si fro aech sceciiafplly difrnfeet"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble_word(scrambled):
    for perm in itertools.permutations(scrambled):
        candidate = ''.join(perm)
        if candidate in word_list:
            return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)