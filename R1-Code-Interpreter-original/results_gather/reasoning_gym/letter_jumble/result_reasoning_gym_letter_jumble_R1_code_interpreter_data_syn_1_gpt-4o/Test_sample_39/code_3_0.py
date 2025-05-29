# The scrambled sentence
scrambled_sentence = "Verne ordRatce s neot eht Year saw frits upblidhes ni teh rbFyuare p"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# A refined list of possible words
possible_words = [
    "Verne", "creator", "as", "note", "the", "Year", "was", "first", "published", "in", "February", "by"
]

# Function to unscramble a word
def unscramble_word(scrambled_word, possible_words):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Check against each possible word
    for word in possible_words:
        if sorted_scrambled == sorted(word.lower()):
            # Preserve the original capitalization
            if scrambled_word.istitle():
                return word.title()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word, possible_words) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)