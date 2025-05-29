import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "forgot", "to", "note", "the", "time", "and", "so", "am", "a", "little", "tale", "yes", "a", "little", "world"
]

# The scrambled sentence
scrambled_sentence = "oogtrf ot note het tmei nda os ma a tetlli tlae Yes a telitl grelowd"

# Function to unscramble a word
def unscramble_word(scrambled_word, dictionary):
    sorted_scrambled = sorted(scrambled_word.lower())
    for word in dictionary:
        if sorted_scrambled == sorted(word):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble the sentence
unscrambled_words = []
for word in scrambled_sentence.split():
    unscrambled_word = unscramble_word(word, dictionary)
    # Preserve capitalization
    if word[0].isupper():
        unscrambled_word = unscrambled_word.capitalize()
    unscrambled_words.append(unscrambled_word)

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)
print(unscrambled_sentence)