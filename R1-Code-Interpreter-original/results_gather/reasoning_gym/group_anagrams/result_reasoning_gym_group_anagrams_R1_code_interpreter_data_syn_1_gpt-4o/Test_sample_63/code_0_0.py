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

words = ["glowered", "reglowed", "overs", "servo", "verso", "needham", "headmen", 
         "trollopian", "pollinator", "soir", "sori", "astatizer", "zaratites", 
         "gallian", "lingala", "triangulated", "adulterating", "autoradiographic", 
         "radioautographic", "requitals", "quartiles"]

result = group_anagrams(words)
print(result)