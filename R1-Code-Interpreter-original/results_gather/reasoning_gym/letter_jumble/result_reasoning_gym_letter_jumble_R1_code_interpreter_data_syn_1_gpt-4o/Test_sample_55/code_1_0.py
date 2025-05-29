import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# Use the nltk words corpus as the dictionary
dictionary = set(words.words())

# Scrambled words
scrambled_words = ["tish", "norginm", "was", "tno", "giragdn", "fo", "gbeilsnss"]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = sorted(word)
        for dict_word in dictionary:
            if sorted(dict_word) == sorted_word:
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

print(final_sentence)