from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        # Sort the word to form the key
        sorted_word = ''.join(sorted(word))
        # Append the word to the corresponding list in the dictionary
        anagrams[sorted_word].append(word)
    # Return the grouped anagrams as a list of lists
    return list(anagrams.values())

words = ["mal", "lam", "nonvisitation", "innovationist", "sept", "pest", "step", "pets", "osmophore", "sophomore", "malleted", "metalled"]
result = group_anagrams(words)
print(result)