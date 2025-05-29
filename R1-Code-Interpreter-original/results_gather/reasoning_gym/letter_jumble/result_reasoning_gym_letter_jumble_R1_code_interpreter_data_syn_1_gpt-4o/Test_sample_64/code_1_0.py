# List of scrambled words
scrambled_words = [
    "attacehd", "lufl", "rPocejt", "Gbtenuerg", "enceisL", "whne", "yuo", 
    "sahre", "ti", "twthoui", "cahgre", "hitw", "rtoehs", "D"
]

# A small set of possible correct words
possible_words = [
    "attached", "full", "Project", "Gutenberg", "License", "when", "you", 
    "share", "it", "without", "charge", "with", "others", "D"
]

# Function to unscramble a word
def unscramble(word, possible_words):
    word_lower = word.lower()
    for possible_word in possible_words:
        if sorted(word_lower) == sorted(possible_word.lower()):
            # Return the word with the original capitalization
            return possible_word.capitalize() if word[0].isupper() else possible_word
    return word  # Return the original word if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word, possible_words) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)