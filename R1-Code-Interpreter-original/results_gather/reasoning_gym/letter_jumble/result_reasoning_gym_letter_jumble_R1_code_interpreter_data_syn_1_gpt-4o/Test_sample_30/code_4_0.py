# Mapping of scrambled words to their correct forms
word_map = {
    "dhsi": "this",
    "Fro": "For",
    "lkie": "like",
    "lla": "all",
    "aewlthy": "wealthy",
    "oflk": "folk",
    "ni": "in",
    "uor": "our"
}

# List of scrambled words
scrambled_words = ["dhsi", "Fro", "lkie", "lla", "aewlthy", "oflk", "ni", "uor"]

# Unscramble each word using the map
unscrambled_sentence = ' '.join(word_map[word] for word in scrambled_words)

# Print the unscrambled sentence
print(unscrambled_sentence)