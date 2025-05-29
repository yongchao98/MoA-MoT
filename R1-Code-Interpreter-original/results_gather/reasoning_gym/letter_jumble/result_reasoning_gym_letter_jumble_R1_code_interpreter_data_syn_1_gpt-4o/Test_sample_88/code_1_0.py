# List of scrambled words
scrambled_words = ["FO", "DMSAGEA", "Extepc", "ofr", "hte", "Rgiht", "fo", "Rapelcement", "ro", "eRfund", "edscriebd", "ni", "apragrhpa", "F", "hte", "Project", "Guetnbgre", "Liratery", "Arevihc"]

# A small dictionary of common words that might match the scrambled words
common_words = {
    "fo": "of",
    "dmsagea": "damages",
    "extepc": "except",
    "ofr": "for",
    "hte": "the",
    "rgiht": "right",
    "rapelcement": "replacement",
    "ro": "or",
    "erfund": "refund",
    "edscriebd": "described",
    "ni": "in",
    "apragrhpa": "paragraph",
    "project": "project",
    "guetnbgre": "gutenberg",
    "liratery": "literary",
    "arevihc": "archive"
}

# Function to unscramble a word using the common words dictionary
def unscramble_word(scrambled):
    scrambled_lower = scrambled.lower()
    return common_words.get(scrambled_lower, scrambled)

# Unscramble each word
unscrambled_sentence = ' '.join(unscramble_word(word) for word in scrambled_words)

print(unscrambled_sentence)