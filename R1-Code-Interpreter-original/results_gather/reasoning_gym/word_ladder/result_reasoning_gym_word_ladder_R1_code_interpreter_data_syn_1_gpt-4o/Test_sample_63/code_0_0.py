import nltk
from nltk.corpus import words

# Ensure the words corpus is downloaded
nltk.download('words')

# List of valid English words
word_list = set(words.words())

# Function to check if a word is valid
def is_valid_word(word):
    return word.lower() in word_list

# Possible transformations
transformations = [
    "SEED", "BEED", "BADE"
]

# Check the validity of each transformation
valid_sequence = all(is_valid_word(word) for word in transformations)

# Print the result
print(valid_sequence)