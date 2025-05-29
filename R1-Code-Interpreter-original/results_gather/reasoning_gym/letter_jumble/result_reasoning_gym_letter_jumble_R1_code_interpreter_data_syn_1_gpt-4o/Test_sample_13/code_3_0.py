# A more comprehensive dictionary of words
dictionary = [
    "understands", "who", "a", "man", "situation", "as", "Smith", "is", "must", "be", "best"
]

# The scrambled sentence
scrambled_sentence = "asdenutrnds ohw a amn stduetai sa mSiht si msut eb tesbe"

# Function to unscramble words
def unscramble(scrambled_word, dictionary):
    sorted_scrambled = sorted(scrambled_word.lower())
    for word in dictionary:
        if sorted_scrambled == sorted(word.lower()):
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble the sentence
unscrambled_words = [
    unscramble(word, dictionary) for word in scrambled_sentence.split()
]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)