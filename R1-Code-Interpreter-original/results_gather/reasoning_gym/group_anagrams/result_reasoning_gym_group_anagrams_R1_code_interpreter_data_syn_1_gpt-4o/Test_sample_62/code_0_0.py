def group_anagrams(words):
    from collections import defaultdict
    
    anagrams = defaultdict(list)
    
    for word in words:
        # Sort the word to form the key
        key = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[key].append(word)
    
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["siblings", "bislings", "spane", "aspen", "snape", "spean", "sceptering", 
         "respecting", "coil", "coli", "clio", "loci", "contractu", "cunctator", 
         "enduringly", "underlying", "spelaean", "seaplane", "electroacoustic", 
         "acoustoelectric", "moistener", "neoterism", "gartering", "regrating"]

result = group_anagrams(words)
print(result)