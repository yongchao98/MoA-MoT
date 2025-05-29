import itertools
import nltk
from nltk.corpus import words

# Download the words corpus if not already downloaded
nltk.download('words')

# List of scrambled words
scrambled_words = [
    "rutenbeGg", "rmaketard", "ubt", "eh", "ahs", "ageder", "ot", 
    "eonadt", "rlyaoties", "erdun", "shit", "pgraaraph", "ot", 
    "hte", "Pcojert", "Gunteberg", "Liretyar", "Aechriv", "tounnFoiad"
]

# Get the list of English words
english_words = set(words.words())

# Function to unscramble a word
def unscramble_word(scrambled):
    scrambled_lower = scrambled.lower()
    for perm in itertools.permutations(scrambled_lower):
        candidate = ''.join(perm)
        if candidate in english_words:
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)