import enchant

# Initialize the English dictionary
d = enchant.Dict("en_US")

# Scrambled words
scrambled_words = ["wtleve", "eth", "uroh"]

# Possible unscrambled words
unscrambled_words = ["twelve", "the", "hour"]

# Check if the unscrambled words are valid
valid_words = [d.check(word) for word in unscrambled_words]

# Print the result
print(all(valid_words))