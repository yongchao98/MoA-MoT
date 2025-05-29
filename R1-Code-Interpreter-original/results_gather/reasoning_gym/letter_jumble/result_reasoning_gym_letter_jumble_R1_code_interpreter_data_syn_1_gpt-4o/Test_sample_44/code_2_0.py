# Predefined dictionary of words
dictionary = {
    "United": ["tinUed", "dniUet"],
    "States": "Sattes",
    "and": "nad",
    "you": "oyu",
    "are": "rae",
    "located": "dotecal",
    "in": "ni",
    "the": "teh",
    "we": "ew",
    "do": "od",
    "not": "ton",
    "claim": "cmali"
}

# List of scrambled words
scrambled_sentence = "tinUed Sattes nad oyu rae dotecal ni teh dniUet States ew od ton cmali"
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word using the dictionary
def unscramble_word(scrambled):
    for word, scrambled_forms in dictionary.items():
        if isinstance(scrambled_forms, list):
            if scrambled.lower() in map(str.lower, scrambled_forms):
                return word
        elif scrambled.lower() == scrambled_forms.lower():
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)