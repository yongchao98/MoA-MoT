# List of scrambled words
scrambled_words = ["nncldiiug", "pagiyn", "yoasltier", "fro", "esu", "fo", "teh", "jPorect"]

# Manually defined dictionary of possible words
possible_words = {
    "nncldiiug": "including",
    "pagiyn": "paying",
    "yoasltier": "royalties",
    "fro": "for",
    "esu": "use",
    "fo": "of",
    "teh": "the",
    "jPorect": "Project"
}

# Unscramble each word using the dictionary
unscrambled_words = [possible_words[word] for word in scrambled_words]

# Reconstruct the sentence
sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(sentence)