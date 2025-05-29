from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to get the key
        key = ''.join(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[key].append(word)
    # Return the values of the dictionary as a list of lists
    return list(anagrams.values())

words = ["calorist", "coralist", "unsoiling", "insouling", "sticklers", "strickles", "trickless"]
result = group_anagrams(words)
print(result)