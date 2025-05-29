# A predefined list of words that includes the necessary words
dictionary = {"of", "the", "people", "Great", "Britain", "who", "still", "now", "had"}

# Scrambled words
scrambled_words = ["fo", "teh", "eeoplp", "fo", "areGt", "Birtain", "hwo", "tlli", "nwo", "ahd"]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for word in dictionary:
        if sorted(word.lower()) == sorted(scrambled.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

print(unscrambled_sentence)