import itertools

# A small dictionary of words for demonstration purposes
dictionary = [
    "electronic", "work", "or", "group", "of", "works", "on", "different", 
    "terms", "than", "are", "set", "forth", "in"
]

# Scrambled words
scrambled_words = [
    "electconir", "wokr", "ro", "proug", "fo", "sorkw", "no", 
    "tifferend", "temrs", "htan", "aer", "est", "rofth", "ni"
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

# Join the unscrambled words into a sentence
unscrambled_sentence_str = ' '.join(unscrambled_sentence)

# Print the unscrambled sentence
print(unscrambled_sentence_str)