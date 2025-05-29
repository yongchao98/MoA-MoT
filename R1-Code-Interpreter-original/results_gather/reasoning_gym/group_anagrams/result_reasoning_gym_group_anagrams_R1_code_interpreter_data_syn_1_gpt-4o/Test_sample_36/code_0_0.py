words = ["nomological", "monological", "predelay", "replayed", "parleyed", "restretch", "stretcher", "unpooled", "unlooped", "gein", "gien", "arginine", "nigerian", "makars", "karmas", "inhales", "leanish", "cenchrus", "crunches"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)