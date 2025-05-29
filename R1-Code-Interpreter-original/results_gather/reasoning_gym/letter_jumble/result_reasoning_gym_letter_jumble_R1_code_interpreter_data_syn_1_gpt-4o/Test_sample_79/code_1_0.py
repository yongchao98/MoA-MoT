# List of scrambled words
scrambled_words = ["teeeohpnl", "I", "was", "wtih", "rD", "Wlkisni", "hA", "the", "dogo"]

# A simple dictionary of possible correct words
# This is a limited set for demonstration purposes
simple_dictionary = {
    "teeeohpnl": "telephone",
    "wtih": "with",
    "rD": "Dr",
    "Wlkisni": "Wilkins",
    "hA": "Ah",
    "dogo": "good"
}

# Function to unscramble a word using the simple dictionary
def unscramble(word):
    return simple_dictionary.get(word, word)

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

print(unscrambled_sentence)