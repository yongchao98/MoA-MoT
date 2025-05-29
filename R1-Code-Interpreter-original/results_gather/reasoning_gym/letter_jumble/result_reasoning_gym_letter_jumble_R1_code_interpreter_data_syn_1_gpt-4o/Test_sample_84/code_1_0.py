# List of scrambled words
scrambled_words = [
    "arYe", "asw", "tirsf", "eshlibpdu", "ni", "eth", "rebuaryF", "p", 
    "tI", "asw", "epsibluhd", "ni", "Fnerca", "het", "tnxe", "year", 
    "thghouAl", "pdlbsihue", "nrued", "eth"
]

# Manually defined dictionary of possible words
possible_words = [
    "Year", "was", "first", "published", "in", "the", "February", "by", 
    "It", "was", "published", "in", "France", "the", "next", "year", 
    "Although", "published", "under", "the"
]

# Function to unscramble a word
def unscramble_word(scrambled, possible_words):
    scrambled_lower = scrambled.lower()
    for word in possible_words:
        if sorted(scrambled_lower) == sorted(word.lower()):
            return word
    return scrambled

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word, possible_words) for word in scrambled_words]

# Reconstruct the sentence
sentence = " ".join(unscrambled_words)

# Print the unscrambled sentence
print(sentence)