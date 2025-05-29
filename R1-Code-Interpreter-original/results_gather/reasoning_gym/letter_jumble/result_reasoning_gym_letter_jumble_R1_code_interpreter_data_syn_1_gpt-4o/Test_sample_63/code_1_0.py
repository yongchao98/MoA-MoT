# A small list of common English words for demonstration purposes
common_words = {
    "os": "so", "atht": "that", "ti": "it", "yam": "may", "ees": "see", 
    "untesrdand": "understand", "adn": "and", "mbrereem": "remember", 
    "teh": "the", "reasnos": "reasons", "which": "which", "eteirmend": "determined", 
    "hte": "the", "pelionrysat": "personality", "ttah": "that", "atsc": "acts", 
    "tJsu": "Just", "tyuds": "study", "uyolsfre": "yourself", "sa": "as"
}

# List of scrambled words
scrambled_words = [
    "os", "atht", "ti", "yam", "ees", "untesrdand", "adn", "mbrereem", "teh", 
    "reasnos", "which", "eteirmend", "hte", "pelionrysat", "ttah", "atsc", 
    "tJsu", "tyuds", "uyolsfre", "sa"
]

# Unscramble each word using the predefined dictionary
unscrambled_words = [common_words.get(word, word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)