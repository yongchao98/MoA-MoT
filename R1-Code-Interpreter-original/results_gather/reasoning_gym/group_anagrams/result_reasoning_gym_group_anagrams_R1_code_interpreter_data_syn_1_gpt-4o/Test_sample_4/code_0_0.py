def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word into a tuple of characters
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[sorted_word].append(word)
    
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["almanac", "mancala", "scaupers", "acupress", "spermalist", "slipstream", 
         "psalmister", "wheeled", "wheedle", "unplied", "unpiled", "praepostor", 
         "pterospora", "coat", "taco"]

result = group_anagrams(words)
print(result)