def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word into a tuple of characters
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[sorted_word].append(word)
    
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["granitite", "iterating", "helium", "humlie", "nizam", "nazim", 
         "striplings", "slipstring", "rearrest", "arrester", "bf", "fb", 
         "tadpolism", "diplomats", "cunan", "canun"]

result = group_anagrams(words)
print(result)