# List of scrambled words
scrambled_words = [
    "nuteGberg", "Ltreairy", "Archeiv", "oouidatnFn", "Pcorejt",
    "uGterbeng", "edpdnes", "npou", "nad", "oancnt", "serivvu",
    "uithowt", "rieewpsdad"
]

# Updated dictionary of possible words
dictionary = [
    "Gutenberg", "Literary", "Archive", "Foundation", "Project",
    "Gutenberg", "depends", "upon", "and", "cannot", "survive",
    "without", "redistribute"
]

# Function to unscramble a word
def unscramble(scrambled, dictionary):
    scrambled_sorted = sorted(scrambled.lower())
    for word in dictionary:
        if scrambled_sorted == sorted(word.lower()):
            return word
    return scrambled

# Debugging: Check the sorted version of "rieewpsdad"
print("Sorted 'rieewpsdad':", sorted("rieewpsdad".lower()))
print("Sorted 'redistribute':", sorted("redistribute".lower()))

# Unscramble each word
unscrambled_sentence = [unscramble(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_sentence)

print(unscrambled_sentence)