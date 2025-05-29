# Expanded dictionary of words relevant to the sentence
expanded_words = {
    "they", "are", "ever", "disregards", "to", "the", "American", 
    "government", "vain", "compilation", "A", "campaign", "against", 
    "the", "annexation", "that", "has"
}

# List of scrambled words
scrambled_words = [
    "tyeh", "are", "erev", "aseridgnsd", "ot", "het", "nearicAm", 
    "nmveonergt", "ivna", "atoilmcspn", "A", "imcpaagn", "stgiaan", 
    "eth", "anaaenxnito", "thatt", "has"
]

# Function to unscramble a word
def unscramble_word(scrambled):
    scrambled_lower = scrambled.lower()
    for word in expanded_words:
        if sorted(word.lower()) == sorted(scrambled_lower):
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Reconstruct the sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)