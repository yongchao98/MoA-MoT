# List of scrambled words
scrambled_words = ["ubt", "fo", "seytyedar", "A", "vulalbae", "enointivn", "ineedd", "dna", "rM", "miSth", "shit", "nnomirg", "wsa", "nto", "nggiard", "fo"]

# Updated dictionary of possible words
common_words = [
    "but", "of", "yesterday", "a", "valuable", "invention", "indeed", "and", 
    "Mr", "Smith", "this", "morning", "was", "not", "dragging", "for"
]

# Function to unscramble a word
def unscramble_word(scrambled, word_list):
    scrambled_sorted = sorted(scrambled.lower())
    for word in word_list:
        if sorted(word.lower()) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = []
for word in scrambled_words:
    unscrambled_word = unscramble_word(word, common_words)
    # Preserve capitalization
    if word.istitle():
        unscrambled_word = unscrambled_word.capitalize()
    unscrambled_sentence.append(unscrambled_word)

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_sentence)

print(unscrambled_sentence)