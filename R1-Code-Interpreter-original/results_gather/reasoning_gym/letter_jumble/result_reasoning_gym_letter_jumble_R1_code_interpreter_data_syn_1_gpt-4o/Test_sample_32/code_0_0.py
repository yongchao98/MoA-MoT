import itertools

# List of scrambled words
scrambled_words = [
    "lulf", "Prcjeot", "Guetnberg", "eicLnse", "tusm", "rppeaa", 
    "yeominrntlp", "wreneveh", "yna", "pocy", "fo", "a", "trojecP", 
    "Gutenbreg", "wkro", "yna", "wkro", "no"
]

# A small dictionary of words for demonstration purposes
# In a real scenario, this would be a comprehensive dictionary
dictionary = [
    "full", "Project", "Gutenberg", "License", "must", "appear", 
    "prominently", "whenever", "any", "copy", "of", "a", "Project", 
    "Gutenberg", "work", "any", "work", "on"
]

# Function to unscramble a word
def unscramble(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled.lower()) == sorted(word.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble(word, dictionary) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)