from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["snicher", "richens", "mathetic", "thematic", "leanish", "inhales", "miasmal", "lamaism", 
         "sencion", "conines", "connies", "romancist", "romantics", "narcotism", "faired", "frieda", 
         "unlighted", "undelight", "elong", "longe"]

result = group_anagrams(words)
print(result)