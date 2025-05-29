from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["radical", "cardial", "evened", "vendee", "reduit", "erudit", "when", "hewn", "allower", "reallow", "apostle", "pelotas", "aseptol"]
result = group_anagrams(words)
print(result)