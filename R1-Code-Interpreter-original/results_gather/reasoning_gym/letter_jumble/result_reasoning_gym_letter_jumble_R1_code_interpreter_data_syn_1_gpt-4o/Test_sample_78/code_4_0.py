# A more comprehensive dictionary of known words
known_words = [
    "go", "to", "note", "the", "time", "and", "so", "am", "a", "little", "tale", "yes", "glowed", "forgot"
]

# Function to unscramble a word
def unscramble_word(scrambled_word, known_words):
    for word in known_words:
        if sorted(scrambled_word.lower()) == sorted(word.lower()):
            return word
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "oogtrf ot note het tmei nda os ma a tetlli tlae Yes a telitl grelowd"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, known_words) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)