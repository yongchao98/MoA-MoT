from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word into a tuple of characters
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[sorted_word].append(word)
    # Return the list of anagram groups
    return list(anagrams.values())

words = ["encrusted", "undersect", "uncrested", "pretaste", "tapester", "prestate", 
         "peartest", "topological", "optological", "chapote", "apothec", "amelcorn", 
         "cornmeal", "ted", "det", "umlauts", "mutuals", "extirpates", "sexpartite"]

result = group_anagrams(words)
print(result)