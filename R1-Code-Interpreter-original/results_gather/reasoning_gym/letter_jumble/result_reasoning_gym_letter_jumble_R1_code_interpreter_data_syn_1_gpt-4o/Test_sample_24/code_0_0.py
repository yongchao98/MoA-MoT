import itertools

# Scrambled words from the sentence
scrambled_words = [
    "ogr", "denato", "ceStion", "reneGal", "Infomration", "Atoub", 
    "erojPct", "eutGnberg", "electronci", "wkros", "orPfessor", 
    "Mhciael", "S", "aHrt", "saw", "eht", "ortginaior", "fo", "eht", "Ptojecr"
]

# A small dictionary of words for demonstration purposes
# In a real scenario, this would be a comprehensive dictionary
dictionary = [
    "org", "donate", "Section", "General", "Information", "About", 
    "Project", "Gutenberg", "electronic", "works", "Professor", 
    "Michael", "S", "Hart", "was", "the", "originator", "of", "the", "Project"
]

# Function to unscramble a word
def unscramble_word(scrambled, dictionary):
    for word in dictionary:
        if sorted(scrambled.lower()) == sorted(word.lower()):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble the sentence
unscrambled_sentence = ' '.join(unscramble_word(word, dictionary) for word in scrambled_words)

print(unscrambled_sentence)