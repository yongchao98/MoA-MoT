def group_anagrams(words):
    anagrams = {}
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Group words by their sorted key
        if key in anagrams:
            anagrams[key].append(word)
        else:
            anagrams[key] = [word]
    # Return the grouped anagrams
    return list(anagrams.values())

words = ["ensuing", "gunnies", "delapse", "elapsed", "sepaled", "iridous", "isuroid", "beatille", "etabelli"]
result = group_anagrams(words)
print(result)