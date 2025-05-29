# List of scrambled words
scrambled_words = ["yma", "seu", "hist", "koeoB", "orf", "elyanr", "any", "ppuoers", "chus", "sa", "oitrcean", "fo", "evrtaiidve", "roswk", "retsorp", "psmoenrracfe", "dan", "erchrase"]

# A more comprehensive dictionary of words
dictionary = [
    "my", "use", "this", "book", "for", "nearly", "any", "purpose", "such", "as", "creation", "of", "derivative", "works", "presort", "performance", "and", "research"
]

# Function to unscramble a word
def unscramble(scrambled, dictionary):
    scrambled_sorted = sorted(scrambled.lower())
    for word in dictionary:
        if scrambled_sorted == sorted(word.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word, dictionary) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = " ".join(unscrambled_words)

print(unscrambled_sentence)