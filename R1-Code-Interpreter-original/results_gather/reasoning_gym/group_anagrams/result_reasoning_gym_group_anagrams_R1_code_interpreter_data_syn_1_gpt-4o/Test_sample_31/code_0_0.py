from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        key = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[key].append(word)
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["epistatic", "pistacite", "commissioned", "decommission", "unportrayed", 
         "unpredatory", "chandelle", "chanelled", "footpad", "padfoot", "noma", 
         "oman", "noam", "innodate", "antinode", "anointed", "distater", 
         "tardiest", "striated", "shapka", "pashka"]

result = group_anagrams(words)
print(result)