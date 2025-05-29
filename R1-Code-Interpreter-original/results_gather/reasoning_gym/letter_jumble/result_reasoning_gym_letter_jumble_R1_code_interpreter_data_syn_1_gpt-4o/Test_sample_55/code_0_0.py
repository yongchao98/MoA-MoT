import itertools

# A small dictionary of words for demonstration purposes
dictionary = ["this", "morning", "was", "not", "grading", "of", "blessing"]

# Scrambled words
scrambled_words = ["tish", "norginm", "was", "tno", "giragdn", "fo", "gbeilsnss"]

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

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

print(final_sentence)