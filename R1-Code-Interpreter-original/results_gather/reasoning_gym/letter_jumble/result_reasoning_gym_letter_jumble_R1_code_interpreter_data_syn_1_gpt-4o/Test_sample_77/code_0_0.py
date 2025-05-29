import itertools

# List of scrambled words
scrambled_words = ["yna", "omney", "pdia", "rof", "a", "wkro", "ro", "a", "retleaemcnp", "pocy", "fi"]

# A small dictionary of possible words (for demonstration purposes)
# In a real scenario, this would be a comprehensive dictionary
dictionary = ["any", "money", "paid", "for", "a", "work", "or", "a", "replacement", "copy", "if"]

# Function to unscramble a word
def unscramble(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled) == sorted(word):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_sentence = [unscramble(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

print(final_sentence)