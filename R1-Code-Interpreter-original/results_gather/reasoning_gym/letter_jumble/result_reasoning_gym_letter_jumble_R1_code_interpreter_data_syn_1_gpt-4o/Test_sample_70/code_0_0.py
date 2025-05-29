import itertools

# A small dictionary of words for demonstration purposes
dictionary = ["United", "States", "without", "permission", "and", "without", "paying"]

# The scrambled sentence
scrambled_sentence = "Utiedn tteaSs wtthoiu erpmisinos dna iwtohut ngiypa"

# Function to unscramble a word
def unscramble_word(scrambled_word, dictionary):
    scrambled_sorted = sorted(scrambled_word.lower())
    for word in dictionary:
        if sorted(word.lower()) == scrambled_sorted:
            return word
    return scrambled_word  # Return the original if no match is found

# Unscramble the sentence
unscrambled_words = []
for word in scrambled_sentence.split():
    unscrambled_word = unscramble_word(word, dictionary)
    unscrambled_words.append(unscrambled_word)

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)