from itertools import permutations

# List of scrambled words
scrambled_words = [
    "nuteGberg", "Ltreairy", "Archeiv", "oouidatnFn", "Pcorejt",
    "uGterbeng", "edpdnes", "npou", "nad", "oancnt", "serivvu",
    "uithowt", "rieewpsdad"
]

# A simple dictionary of possible words
# In a real scenario, this would be a comprehensive dictionary
dictionary = [
    "Gutenberg", "Literary", "Archive", "Foundation", "Project",
    "Gutenberg", "depends", "upon", "and", "cannot", "survive",
    "without", "redistribute"
]

# Function to unscramble a word
def unscramble(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled.lower()) == sorted(word.lower()):
            return word
    return scrambled

# Unscramble each word
unscrambled_sentence = [unscramble(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_sentence)

print(unscrambled_sentence)