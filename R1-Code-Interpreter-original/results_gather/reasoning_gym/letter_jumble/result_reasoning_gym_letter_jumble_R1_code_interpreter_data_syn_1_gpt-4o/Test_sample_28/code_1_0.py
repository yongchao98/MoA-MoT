# Predefined list of possible correct words
possible_words = [
    "broke", "in", "Mr", "Smith", "turning", "to", "a", "second", 
    "scientist", "Optical", "science", "defective"
]

# List of scrambled words
scrambled_words = ["borke", "ni", "rM", "mSith", "ntruing", "ot", "a", "sdcone", "ncseitist", "tpOilac", "cniesce", "defcetive"]

# Function to unscramble a word
def unscramble_word(scrambled, possible_words):
    scrambled_sorted = sorted(scrambled.lower())
    for word in possible_words:
        if sorted(word.lower()) == scrambled_sorted:
            # Preserve the original capitalization
            if scrambled[0].isupper():
                return word.capitalize()
            else:
                return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word, possible_words) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)