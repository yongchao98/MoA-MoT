# Predefined list of valid English words
valid_words = {"twelve", "the", "hour"}

# Scrambled words and their possible unscrambled forms
scrambled_to_unscrambled = {
    "wtleve": "twelve",
    "eth": "the",
    "uroh": "hour"
}

# Check if the unscrambled words are valid
unscrambled_words = [scrambled_to_unscrambled[word] for word in scrambled_to_unscrambled]
are_valid_words = all(word in valid_words for word in unscrambled_words)

# Print the result
print(are_valid_words)