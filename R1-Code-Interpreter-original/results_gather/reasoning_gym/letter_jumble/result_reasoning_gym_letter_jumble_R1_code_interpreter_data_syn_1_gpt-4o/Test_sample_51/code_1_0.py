import itertools

# List of scrambled words
scrambled_words = [
    "impncloete", "enrctaiuca", "ro", "cpturor", "tada", "nrcstiaripton",
    "rerrso", "a", "horiypgct", "ro", "tehro", "autellcletin", "reorppty",
    "tefmeninirgn", "a", "dteecivef", "ro", "degamda"
]

# Updated dictionary of possible words
dictionary = [
    "incomplete", "inaccurate", "or", "corrupt", "data", "transcription",
    "errors", "a", "copyright", "or", "other", "intellectual", "property",
    "infringement", "a", "defective", "or", "damaged"
]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    sorted_scrambled = sorted(scrambled)
    for word in dictionary:
        if sorted(word) == sorted_scrambled:
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

print(unscrambled_sentence)