# List of scrambled words
scrambled_words = ["lovens", "Addrgnsise", "oen", "fo", "teesh", "uathors", "hwo", "aws", "wgitina", "sih", "trun", "Capitla"]

# Manually defined dictionary of possible words
possible_words = {
    "lovens": "novels",
    "Addrgnsise": "Sandgrisen",  # Assuming this is a proper noun
    "oen": "one",
    "fo": "of",
    "teesh": "these",
    "uathors": "authors",
    "hwo": "who",
    "aws": "was",
    "wgitina": "waiting",
    "sih": "his",
    "trun": "turn",
    "Capitla": "Capital"  # Assuming this is a proper noun
}

# Unscramble each word using the predefined dictionary
unscrambled_sentence = ' '.join(possible_words.get(word, word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)