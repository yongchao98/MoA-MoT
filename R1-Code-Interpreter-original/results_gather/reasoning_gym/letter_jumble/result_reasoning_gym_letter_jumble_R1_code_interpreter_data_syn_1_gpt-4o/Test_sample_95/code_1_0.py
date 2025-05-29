import itertools

# A small dictionary of words for demonstration purposes
dictionary = {"the", "first", "thing", "that"}

def unscramble_word(scrambled_word, dictionary):
    # Generate all permutations of the scrambled word
    permutations = set(itertools.permutations(scrambled_word.lower()))
    # Check each permutation to see if it's a valid word
    for perm in permutations:
        candidate = ''.join(perm)
        if candidate in dictionary:
            # Return the word with the original capitalization
            if scrambled_word[0].isupper():
                return candidate.capitalize()
            else:
                return candidate
    return scrambled_word  # Return the original if no match is found

# The scrambled sentence
scrambled_sentence = "heT isrft nhigt hatt"

# Unscramble each word in the sentence
unscrambled_words = [unscramble_word(word, dictionary) for word in scrambled_sentence.split()]

# Join the unscrambled words into a sentence
unscrambled_sentence = ' '.join(unscrambled_words)

print(unscrambled_sentence)