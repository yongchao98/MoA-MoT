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

words = ["sycones", "coyness", "cogniser", "coresign", "coreigns", "glumose", "lugsome", 
         "stanitsa", "satanist", "mistrace", "camister", "camerist", "poietic", "epiotic", 
         "rouse", "euros", "roues", "woodhack", "hackwood"]

result = group_anagrams(words)
print(result)