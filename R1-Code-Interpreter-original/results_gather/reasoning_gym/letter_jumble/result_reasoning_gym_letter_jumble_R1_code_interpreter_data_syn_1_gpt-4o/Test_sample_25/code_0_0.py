import itertools

# Define a function to check if a word is valid
def is_valid_word(word, dictionary):
    return word in dictionary

# A simple dictionary of words
dictionary = {"would", "but", "put", "picture", "to"}

# Scrambled words
scrambled_words = ["wuold", "btu", "piceutr", "ot"]

# Unscramble the words
unscrambled_words = []
for scrambled in scrambled_words:
    for perm in itertools.permutations(scrambled):
        possible_word = ''.join(perm)
        if is_valid_word(possible_word, dictionary):
            unscrambled_words.append(possible_word)
            break

# Join the unscrambled words to form a sentence
sentence = ' '.join(unscrambled_words)
print(sentence)