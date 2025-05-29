def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word to form the key
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[sorted_word].append(word)
    
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["oliveness", "novelises", "reeler", "rereel", "anjan", "jnana", "suitings", "tissuing"]
result = group_anagrams(words)
print(result)