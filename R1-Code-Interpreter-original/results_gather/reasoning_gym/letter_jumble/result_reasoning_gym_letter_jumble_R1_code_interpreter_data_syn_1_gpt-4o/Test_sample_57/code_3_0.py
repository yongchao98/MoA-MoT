from itertools import permutations

# A comprehensive dictionary of words
dictionary = ["an", "invention", "needing", "capital", "again"]

def unscramble_word(scrambled, dictionary):
    for perm in permutations(scrambled):
        candidate = ''.join(perm)
        if candidate in dictionary:
            return candidate
    return scrambled  # Return the original if no match is found

# The scrambled sentence
scrambled_sentence = "na troenivn egendni aatipcl nigaa"
scrambled_words = scrambled_sentence.split()

# Unscramble each word
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_words]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)