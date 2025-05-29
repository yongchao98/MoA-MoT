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

words = ["pantod", "dopant", "helistop", "hoplites", "isopleth", "unteem", "umteen", 
         "nietzsche", "zechstein", "smuttier", "trisetum", "gelosin", "longies", "lignose"]

result = group_anagrams(words)
print(result)