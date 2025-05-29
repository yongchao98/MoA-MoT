# A small set of common English words for demonstration purposes
common_words = {
    "of", "interest", "to", "day", "indeed", "we", "have", "answered", 
    "coyote", "the", "elements", "olympus", "are", "definitely"
}

# List of scrambled words
scrambled_words = [
    "fo", "sntreeit", "ot", "ady", "Indeed", "ew", "ahve", "wnsaered", 
    "Coyoel", "ehT", "lenmeets", "fo", "Opsmluy", "era", "tleidinivefy"
]

# Function to unscramble a word
def unscramble(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = ''.join(sorted(scrambled_word.lower()))
    # Check against the common words
    for word in common_words:
        if sorted_scrambled == ''.join(sorted(word)):
            # Preserve capitalization
            if scrambled_word.istitle():
                return word.title()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)