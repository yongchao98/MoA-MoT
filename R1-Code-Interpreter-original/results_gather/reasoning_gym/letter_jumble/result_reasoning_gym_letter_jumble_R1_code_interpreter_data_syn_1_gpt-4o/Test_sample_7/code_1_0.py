# A small dictionary of common words
common_words = {
    "sa": "as",
    "amn": "man",
    "himself": "himself",
    "rtecleicity": "electricity",
    "aws": "was",
    "knwon": "known",
    "reays": "years",
    "aog": "ago",
    "nad": "and"
}

# List of scrambled words
scrambled_words = ["sa", "amn", "himself", "rtecleicity", "aws", "knwon", "reays", "aog", "nad"]

# Unscramble each word using the dictionary
unscrambled_sentence = ' '.join(common_words[word] for word in scrambled_words)

print(unscrambled_sentence)