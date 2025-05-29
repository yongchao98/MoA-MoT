def group_anagrams(words):
    anagrams = {}
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Add the word to the corresponding anagram group
        if key in anagrams:
            anagrams[key].append(word)
        else:
            anagrams[key] = [word]
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["bouffon", "buffoon", "balanites", "basaltine", "stainable", "pactolian", "placation", "pactional"]
result = group_anagrams(words)
print(result)