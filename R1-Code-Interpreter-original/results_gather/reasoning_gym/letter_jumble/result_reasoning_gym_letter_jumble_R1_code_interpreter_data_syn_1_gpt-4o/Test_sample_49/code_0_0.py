import itertools

# List of scrambled words
scrambled_words = [
    "insrgi", "cseldrea", "het", "nma", "ot", "eb", "adde", "deaD", 
    "lemsaicx", "reevy", "eon", "eteprns", "sYe", "sawnesr", "rD", 
    "iisklnW", "aedd", "And", "how"
]

# A simple dictionary of words for demonstration purposes
# In a real scenario, this should be a comprehensive dictionary
dictionary = [
    "rising", "declares", "the", "man", "to", "be", "dead", "Dead", 
    "exclaims", "every", "one", "presents", "Yes", "answers", "Dr", 
    "Wilkins", "dead", "And", "how"
]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = ''.join(sorted(word.lower()))
        for dict_word in dictionary:
            if sorted_word == ''.join(sorted(dict_word.lower())):
                # Preserve capitalization
                if word[0].isupper():
                    unscrambled_sentence.append(dict_word.capitalize())
                else:
                    unscrambled_sentence.append(dict_word)
                break
    return ' '.join(unscrambled_sentence)

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)
print(unscrambled_sentence)