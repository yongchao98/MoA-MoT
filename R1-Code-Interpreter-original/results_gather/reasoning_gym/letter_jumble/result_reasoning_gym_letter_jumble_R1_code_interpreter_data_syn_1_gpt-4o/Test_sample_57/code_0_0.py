import itertools

# A simple dictionary of words for demonstration purposes
dictionary = ["an", "invention", "needing", "capital", "again"]

def unscramble_word(scrambled, dictionary):
    scrambled_sorted = sorted(scrambled)
    for word in dictionary:
        if sorted(word) == scrambled_sorted:
            return word
    return scrambled  # Return the original if no match is found

# The scrambled sentence
scrambled_sentence = "na troenivn egendni aatipcl nigaa"
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)