# List of scrambled words
scrambled_words = ["IIPLMED", "INCLNDIUG", "TUB", "NTO", "LIMIDET", "OT", "AWRRANTIES", "FO", "NERCHAMTABILTIY", "RO", "FITENSS"]

# Predefined list of possible correct words
possible_words = [
    "implied", "including", "but", "not", "limited", "to", "warranties", "of", "merchantability", "or", "fitness"
]

# Function to unscramble a word
def unscramble(scrambled_word, possible_words):
    # Sort the letters of the scrambled word
    sorted_scrambled = sorted(scrambled_word.lower())
    # Find the matching word by comparing sorted letters
    for word in possible_words:
        if sorted_scrambled == sorted(word):
            return word
    return scrambled_word.lower()  # Return the original if no match is found

# Unscramble each word in the list
unscrambled_sentence = ' '.join(unscramble(word, possible_words) for word in scrambled_words)

# Capitalize the first letter of the sentence
unscrambled_sentence = unscrambled_sentence.capitalize()

print(unscrambled_sentence)