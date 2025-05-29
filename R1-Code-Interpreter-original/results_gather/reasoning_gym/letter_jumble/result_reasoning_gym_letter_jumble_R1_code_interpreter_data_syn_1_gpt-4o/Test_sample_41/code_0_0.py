import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "IMPLIED", "INCLUDING", "BUT", "NOT", "LIMITED", "TO", "WARRANTIES", 
    "OF", "MERCHANTABILITY", "OR", "FITNESS"
]

# Scrambled words
scrambled_words = [
    "IIPLMED", "INCLNDIUG", "TUB", "NTO", "LIMIDET", "OT", 
    "AWRRANTIES", "FO", "NERCHAMTABILTIY", "RO", "FITENSS"
]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = sorted(word)
        for dict_word in dictionary:
            if sorted(dict_word) == sorted_word:
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)
print(" ".join(unscrambled_sentence))