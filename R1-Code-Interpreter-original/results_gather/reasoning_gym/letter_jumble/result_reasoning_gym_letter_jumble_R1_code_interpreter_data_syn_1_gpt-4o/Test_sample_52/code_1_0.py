# List of scrambled words
scrambled_words = ["fo", "tehes", "aosrbb", "dan", "oscnende", "the", "invgil", "crofe", "odninaetc", "ni", "eth", "sun"]

# Predefined list of possible correct words
possible_words = ["of", "these", "absorb", "and", "condense", "the", "living", "force", "contained", "in", "the", "sun"]

# Function to unscramble a word
def unscramble(scrambled_word, possible_words):
    for word in possible_words:
        if sorted(scrambled_word) == sorted(word):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word, possible_words) for word in scrambled_words)

print(unscrambled_sentence)