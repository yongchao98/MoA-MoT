# Direct mapping of scrambled words to unscrambled words
scrambled_to_unscrambled = {
    "ym": "my",
    "omclrucatau": "automobile",
    "wrsko": "works",
    "ta": "at",
    "aNagria": "Niagara",
    "Fro": "For",
    "rM": "Mr",
    "htimS": "Smith",
    "ahs": "has",
    "boiatend": "obtained",
    "a": "a",
    "lease": "lease",
    "fo": "of",
    "hte": "the",
    "arget": "great",
    "sflla": "falls",
    "fo": "of"
}

# Unscramble each word using the direct mapping
unscrambled_words = [scrambled_to_unscrambled[word] for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)