import itertools

# A small dictionary of words for demonstration purposes
dictionary = set([
    "horses", "were", "the", "only", "means", "of", "conveyance", "think", "railroads"
])

# The scrambled sentence
scrambled_sentence = "rohses rewe eht oyln menas fo conveyenca Thikn fo hte rialroads fo teh"

# Function to unscramble a word
def unscramble_word(scrambled_word, dictionary):
    for word in dictionary:
        if sorted(word.lower()) == sorted(scrambled_word.lower()):
            # Preserve the original capitalization
            if scrambled_word.istitle():
                return word.capitalize()
            else:
                return word
    return scrambled_word  # Return the original if no match is found

# Split the sentence into words
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)