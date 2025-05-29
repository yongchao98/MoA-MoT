# Updated dictionary with correct words
dictionary = [
    "gutenberg", "rog", "license", "section", "general", "terms", "of", "use", "and"
]

# The scrambled words
scrambled_words = ["nutegberg", "rog", "license", "Sentioc", "rGnleae", "sTrme", "fo", "Ues", "dna"]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    scrambled_lower = scrambled.lower()
    for word in dictionary:
        if sorted(scrambled_lower) == sorted(word):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = []
for word in scrambled_words:
    unscrambled_word = unscramble_word(word, dictionary)
    # Preserve capitalization
    if word[0].isupper():
        unscrambled_word = unscrambled_word.capitalize()
    unscrambled_sentence.append(unscrambled_word)

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_sentence)
print(unscrambled_sentence)