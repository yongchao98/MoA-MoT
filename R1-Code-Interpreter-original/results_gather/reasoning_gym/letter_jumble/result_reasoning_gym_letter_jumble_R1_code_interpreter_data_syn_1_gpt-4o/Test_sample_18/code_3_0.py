# A small list of common English words relevant to the context
common_words = [
    "may", "use", "this", "Book", "for", "nearly", "any", "purpose", 
    "such", "as", "creation", "of", "derivative", "works", "report", 
    "performance", "and", "research"
]

# Create a dictionary mapping sorted character sequences to words
word_dict = {''.join(sorted(word)): word for word in common_words}

# List of scrambled words
scrambled_words = [
    "yma", "seu", "hist", "koeoB", "orf", "elyanr", "any", "ppuoers", 
    "chus", "sa", "oitrcean", "fo", "evrtaiidve", "roswk", "retsorp", 
    "psmoenrracfe", "dan", "erchrase"
]

# Function to unscramble a word using the dictionary
def unscramble(scrambled_word):
    sorted_word = ''.join(sorted(scrambled_word))
    return word_dict.get(sorted_word, scrambled_word)

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)