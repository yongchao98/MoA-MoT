import itertools

# Scrambled words from the sentence
scrambled_words = [
    "het", "hlnetg", "fo", "the", "aetrsrtrile", "aeyr", "lJues", "reVne", 
    "ndE", "fo", "tcojePr", "nGbreuegt", "s", "nI", "het", "eYar", "yb", "Jlseu"
]

# A simple dictionary of words (for demonstration purposes, a real dictionary would be much larger)
dictionary = [
    "the", "length", "of", "the", "terrestrial", "year", "Jules", "Verne", 
    "End", "of", "Project", "Gutenberg", "is", "In", "the", "Year", "by", "Jules"
]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        # Find a match in the dictionary
        match = next((w for w in dictionary if sorted(w.lower()) == sorted(word.lower())), word)
        unscrambled_sentence.append(match)
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

# Print the final unscrambled sentence
print(final_sentence)