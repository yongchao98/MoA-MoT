import itertools

# List of scrambled words
scrambled_words = ["hcolera", "teh", "yewlol", "feevr", "ttah", "eht", "rllbigeeents", "migth", "spdear"]

# A simple dictionary of words for demonstration purposes
# In a real scenario, this would be a comprehensive dictionary
dictionary = {"cholera", "the", "yellow", "fever", "that", "the", "belligerents", "might", "spread"}

# Function to unscramble a word
def unscramble_word(scrambled_word, dictionary):
    for perm in itertools.permutations(scrambled_word):
        candidate = ''.join(perm)
        if candidate in dictionary:
            return candidate
    return scrambled_word  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)