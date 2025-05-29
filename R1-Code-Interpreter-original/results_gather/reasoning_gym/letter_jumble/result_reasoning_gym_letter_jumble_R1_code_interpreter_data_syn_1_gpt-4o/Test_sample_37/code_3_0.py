# A more comprehensive dictionary of words
dictionary = {"definitely", "settled", "That", "great"}

def unscramble_word(scrambled_word):
    # Sort the letters of the scrambled word
    sorted_scrambled = ''.join(sorted(scrambled_word.lower()))
    # Check against each word in the dictionary
    for word in dictionary:
        if sorted_scrambled == ''.join(sorted(word.lower())):
            # Preserve the original capitalization
            return word if scrambled_word[0].islower() else word.capitalize()
    return scrambled_word

# Scrambled sentence
scrambled_sentence = ["definytiveli", "settdel", "Ttah", "greta"]

# Unscramble each word
unscrambled_sentence = [unscramble_word(word) for word in scrambled_sentence]

# Join the unscrambled words into a sentence
result_sentence = ' '.join(unscrambled_sentence)

print(result_sentence)