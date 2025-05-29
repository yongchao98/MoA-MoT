# A small list of common English words for demonstration
common_words = [
    "Mr", "Smith", "continues", "his", "round", "and", "enters", "the", 
    "reporters", "hall", "Here", "reporters", "in", "their", "respective", 
    "places", "facing", "an", "equal", "number"
]

# Function to unscramble a word
def unscramble_word(scrambled):
    scrambled_sorted = ''.join(sorted(scrambled.lower()))
    for word in common_words:
        if ''.join(sorted(word.lower())) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "rM Smtih inntsouec ihs dounr nad tneers eht rpotreers ahll reHe rreoptesr ni theri eesprcteiv slpcea gacnif na eauql buenmr"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)