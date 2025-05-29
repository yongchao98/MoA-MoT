# Predefined dictionary of words
dictionary = {
    "United": "tinUed",
    "States": "Sattes",
    "and": "nad",
    "you": "oyu",
    "are": "rae",
    "located": "dotecal",
    "in": "ni",
    "the": "teh",
    "do": "od",
    "not": "ton",
    "claim": "cmali"
}

# List of scrambled words
scrambled_sentence = "tinUed Sattes nad oyu rae dotecal ni teh dniUet States ew od ton cmali"
scrambled_words = scrambled_sentence.split()

# Function to unscramble a word using the dictionary
def unscramble_word(scrambled):
    for word, scrambled_form in dictionary.items():
        if scrambled.lower() == scrambled_form.lower():
            return word
    return scrambled  # Return the original if no match is found

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)