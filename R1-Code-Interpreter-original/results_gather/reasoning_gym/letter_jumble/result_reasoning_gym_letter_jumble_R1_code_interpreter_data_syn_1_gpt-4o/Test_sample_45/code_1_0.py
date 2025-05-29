# A small sample list of common English words for demonstration purposes
word_list = [
    "sorts", "of", "dishes", "as", "a", "varied", "assortment", "is", "always", "in", "readiness"
]

def unscramble_sentence(scrambled_sentence):
    scrambled_words = scrambled_sentence.split()
    unscrambled_words = []

    for scrambled_word in scrambled_words:
        # Sort the letters of the scrambled word
        sorted_scrambled = ''.join(sorted(scrambled_word))
        
        # Find the correct word by checking anagrams
        for word in word_list:
            if ''.join(sorted(word)) == sorted_scrambled:
                unscrambled_words.append(word)
                break

    return ' '.join(unscrambled_words)

scrambled_sentence = "osrts fo deshsi sa a idrvea astorsnemt si aswyal ni rsaeinesd"
unscrambled_sentence = unscramble_sentence(scrambled_sentence)
print(unscrambled_sentence)