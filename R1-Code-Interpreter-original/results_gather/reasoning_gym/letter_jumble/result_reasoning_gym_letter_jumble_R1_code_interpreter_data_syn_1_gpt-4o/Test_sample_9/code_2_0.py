# Scrambled words
scrambled_words = [
    "fo", "sntreeit", "ot", "ady", "Indeed", "ew", "ahve", "wnsaered", 
    "Coyoel", "ehT", "lenmeets", "fo", "Opsmluy", "era", "tleidinivefy"
]

# Updated dictionary with correct words
dictionary = [
    "of", "interest", "to", "day", "Indeed", "we", "have", "answered", 
    "Coyote", "The", "elements", "of", "Plumose", "are", "definitely"
]

# Function to unscramble words
def unscramble(scrambled, dictionary):
    unscrambled_sentence = []
    for word in scrambled:
        sorted_word = sorted(word.lower())
        for dict_word in dictionary:
            if sorted(dict_word.lower()) == sorted_word:
                unscrambled_sentence.append(dict_word)
                break
    return unscrambled_sentence

# Unscramble the sentence
unscrambled_sentence = unscramble(scrambled_words, dictionary)

# Join the unscrambled words into a sentence
final_sentence = ' '.join(unscrambled_sentence)

# Print the final unscrambled sentence
print(final_sentence)