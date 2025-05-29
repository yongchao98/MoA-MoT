from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to get the key
        sorted_word = ''.join(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[sorted_word].append(word)
    # Return the grouped anagrams
    return list(anagrams.values())

words = ["eternities", "entireties", "ophis", "hopis", "conformer", "reconform", "mf", "fm"]
result = group_anagrams(words)
print(result)