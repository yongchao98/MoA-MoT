import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "novels", "Gardenside", "one", "of", "these", "authors", "who", "was", "waiting", "his", "turn", "Capital"
]

# Scrambled words
scrambled_words = ["lovens", "Addrgnsise", "oen", "fo", "teesh", "uathors", "hwo", "aws", "wgitina", "sih", "trun", "Capitla"]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = sorted(word.lower())
        for dict_word in dictionary:
            if sorted_word == sorted(dict_word.lower()):
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

print(final_sentence)