import itertools

# List of scrambled words
scrambled_words = [
    "na", "eontnsixe", "craih", "nheT", "nouchigt", "a", "kbon", "eh", 
    "hbtalsiseed", "tomumnicacion", "whti", "teh", "Clnrtae", "Concert", 
    "Hall", "nhcwee", "rou", "grtaeest", "sned", "tuo"
]

# Manually defined dictionary of possible words
possible_words = [
    "an", "extension", "chair", "then", "touching", "a", "book", "he", 
    "established", "communication", "with", "the", "Central", "Concert", 
    "Hall", "when", "our", "greatest", "sends", "out"
]

# Function to unscramble a word
def unscramble(scrambled_word, possible_words):
    for word in possible_words:
        if sorted(scrambled_word.lower()) == sorted(word.lower()):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word, possible_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

print(sentence)