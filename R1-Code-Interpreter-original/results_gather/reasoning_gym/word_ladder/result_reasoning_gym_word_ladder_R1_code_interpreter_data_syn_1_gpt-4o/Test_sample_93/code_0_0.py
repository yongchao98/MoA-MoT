import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of words to check
word_list = ["GIVE", "LIVE", "LOVE", "LOPE", "VEEP"]

# Check if all words are valid
valid_words = set(words.words())
valid_sequence = all(word.lower() in valid_words for word in word_list)

print(valid_sequence)