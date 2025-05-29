# Let's verify if the words are valid English words
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of words to verify
word_list = ["SNOG", "SONG", "SUNG", "SUNK", "SUKS", "SUQS"]

# Check if each word is in the English words corpus
valid_words = [word for word in word_list if word.lower() in words.words()]

# Print the valid words
print(valid_words)