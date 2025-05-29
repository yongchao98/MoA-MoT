# The scrambled sentence
scrambled_sentence = "iMehalc S arHt wsa eht iigortaron fo eht jrcPeto Gutnrbege ecntopc fo a iralbry fo eneriloctc skrow"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# List of potential words based on the context
potential_words = [
    "Michael", "S", "Hart", "was", "the", "originator", "of", "the", 
    "Project", "Gutenberg", "concept", "of", "a", "library", "of", 
    "electronic", "works"
]

# Function to unscramble a single word
def unscramble_word(scrambled_word, potential_words):
    # Sort the characters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check each potential word
    for word in potential_words:
        if sorted_scrambled == sorted(word.lower()):
            return word
    # If no match is found, return the original scrambled word
    return scrambled_word

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word, potential_words) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)