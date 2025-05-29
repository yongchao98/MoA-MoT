words = ["omnitenent", "intonement", "rightable", "bergalith", "missary", "syriasm"]

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

# Group the anagrams
result = group_anagrams(words)
print(result)