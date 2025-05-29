# A small list of common English words relevant to the context
common_words = {
    "may", "use", "this", "Book", "for", "nearly", "any", "purpose", 
    "such", "as", "creation", "of", "derivative", "works", "report", 
    "performance", "and", "research"
}

# List of scrambled words
scrambled_words = [
    "yma", "seu", "hist", "koeoB", "orf", "elyanr", "any", "ppuoers", 
    "chus", "sa", "oitrcean", "fo", "evrtaiidve", "roswk", "retsorp", 
    "psmoenrracfe", "dan", "erchrase"
]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in common_words:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)