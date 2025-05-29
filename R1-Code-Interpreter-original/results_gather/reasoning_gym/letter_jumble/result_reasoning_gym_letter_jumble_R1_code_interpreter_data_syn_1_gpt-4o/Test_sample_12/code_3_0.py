# List of scrambled words
scrambled_words = ["dear", "rchAdbali", "keep", "no", "nSice", "yrstyedae", "atnhks", "ot", "oyu", "rheet", "si", "a", "gain", "fo", "sucrbbisesr", "rM"]

# List of common English words
common_words = [
    "dear", "charitable", "keep", "on", "since", "yesterday", "thanks", "to", "you", "there", "is", "a", "gain", "of", "subscribers", "mr"
]

# Function to find unscrambled word
def unscramble(scrambled_word):
    scrambled_word_lower = scrambled_word.lower()
    for word in common_words:
        if sorted(scrambled_word_lower) == sorted(word.lower()):
            # Preserve the original capitalization
            if scrambled_word[0].isupper():
                return word.capitalize()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word) for word in scrambled_words]

# Reconstruct the sentence with the original capitalization and punctuation
unscrambled_sentence = " ".join(unscrambled_words)
print(unscrambled_sentence)