# Predefined list of words to match against
word_list = ["we", "shall", "be", "able", "to", "produce"]

# Function to unscramble a word
def unscramble(scrambled_word):
    for word in word_list:
        if sorted(word) == sorted(scrambled_word):
            return word
    return scrambled_word  # Return the original if no match is found

# Scrambled sentence
scrambled_sentence = "ew hsall eb ebla ot puodrce"

# Unscramble each word in the sentence
unscrambled_words = [unscramble(word) for word in scrambled_sentence.split()]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)