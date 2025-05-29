from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        sorted_word = tuple(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[sorted_word].append(word)
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["theine", "ethine", "generated", "greatened", "unoffensiveness", "noneffusiveness", "fluorated", "outflared", "ramona", "oarman", "romana"]
result = group_anagrams(words)
print(result)