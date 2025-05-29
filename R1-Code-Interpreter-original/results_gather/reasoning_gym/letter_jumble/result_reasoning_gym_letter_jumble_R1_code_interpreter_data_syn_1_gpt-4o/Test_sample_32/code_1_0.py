# List of scrambled words
scrambled_words = [
    "lulf", "Prcjeot", "Guetnberg", "eicLnse", "tusm", "rppeaa", 
    "yeominrntlp", "wreneveh", "yna", "pocy", "fo", "a", 
    "trojecP", "Gutenbreg", "wkro", "yna", "wkro", "no"
]

# Predefined list of potential words
potential_words = [
    "full", "Project", "Gutenberg", "License", "must", "appear", 
    "prominently", "whenever", "any", "copy", "of", "a", 
    "Project", "Gutenberg", "work", "any", "work", "on"
]

# Function to find unscrambled word
def unscramble(scrambled, potential_words):
    scrambled_lower = scrambled.lower()
    for word in potential_words:
        if sorted(word.lower()) == sorted(scrambled_lower):
            # Preserve capitalization
            if scrambled[0].isupper():
                return word.capitalize()
            else:
                return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble(word, potential_words) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)