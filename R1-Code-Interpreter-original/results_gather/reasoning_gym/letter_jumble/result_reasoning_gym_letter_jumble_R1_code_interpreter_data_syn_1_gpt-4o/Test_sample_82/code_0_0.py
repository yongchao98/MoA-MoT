import itertools

# A simple dictionary of words for demonstration purposes
dictionary = {"of", "wood", "They", "are", "now", "playing", "what", "won", "electricity", "No"}

# The scrambled sentence
scrambled_sentence = "fo wodo eThy era onw ngylaipp hate nwo eeritciyctl oN"

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word
def unscramble_word(scrambled_word, dictionary):
    for word in dictionary:
        if sorted(word.lower()) == sorted(scrambled_word.lower()):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)