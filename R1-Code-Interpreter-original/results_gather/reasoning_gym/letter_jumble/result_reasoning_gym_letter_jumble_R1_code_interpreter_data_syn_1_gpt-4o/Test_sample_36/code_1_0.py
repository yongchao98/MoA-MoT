import itertools

# Scrambled words
scrambled_words = [
    "ePctrjo", "gbeunteGr", "Lieratry", "hAvcrie", "ndotnaFiou", 
    "heT", "Pjtoecr", "rGtegenbu", "rLieyrta", "Acrhive", 
    "Fndtoiuaon", "si", "a", "nno", "ritofp", "c", "eicuotnlaad"
]

# Updated dictionary with the missing word
dictionary = [
    "Project", "Gutenberg", "Literary", "Archive", "Foundation", 
    "The", "is", "a", "non", "profit", "educational"
]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled.lower()) == sorted(word.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)