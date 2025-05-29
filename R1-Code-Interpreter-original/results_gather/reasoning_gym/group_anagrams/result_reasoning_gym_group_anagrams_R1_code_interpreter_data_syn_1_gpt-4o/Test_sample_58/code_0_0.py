words = ["rondelles", "norselled", "criminate", "carminite", "metrician", "antimeric", "flook", "kloof", "stamnoi", "manitos", "incisure", "sciurine"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[key].append(word)
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

# Group the anagrams and print the result
result = group_anagrams(words)
print(result)