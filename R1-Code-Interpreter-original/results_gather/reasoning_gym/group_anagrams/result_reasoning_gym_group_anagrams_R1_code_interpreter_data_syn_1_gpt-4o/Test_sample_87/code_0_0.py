from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        key = ''.join(sorted(word))
        # Append the word to the corresponding anagram group
        anagrams[key].append(word)
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["binder", "rebind", "brined", "inbred", "saxonize", "oxazines", "snop", "pons", "wham", "hawm", "misacts", "mastics", "reparel", "pearler"]
result = group_anagrams(words)
print(result)