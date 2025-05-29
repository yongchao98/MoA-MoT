def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word to form the key
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[sorted_word].append(word)
    
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["navet", "nevat", "venta", "runkly", "knurly", "roqueting", "requoting", "archine", "chainer", "scrime", "crimes"]
result = group_anagrams(words)
print(result)