# List of scrambled words
scrambled_words = ["Faithnubr", "saw", "fdfconei", "adn", "ldia", "ni", "a"]

# Manually defined dictionary of possible words
possible_words = {
    "Faithnubr": "Faithburn",
    "saw": "was",
    "fdfconei": "confided",
    "adn": "and",
    "ldia": "laid",
    "ni": "in",
    "a": "a"
}

# Unscramble each word using the predefined dictionary
unscrambled_sentence = ' '.join(possible_words[word] for word in scrambled_words)

print(unscrambled_sentence)