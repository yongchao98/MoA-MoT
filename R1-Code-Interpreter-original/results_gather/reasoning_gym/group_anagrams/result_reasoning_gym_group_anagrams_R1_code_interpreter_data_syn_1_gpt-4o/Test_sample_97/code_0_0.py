words = ["cynosural", "consulary", "maral", "malar", "alarm", "windless", "wildness", "swindles", "urocele", "urceole", "tawsing", "wasting", "mucoids", "muscoid", "drafters", "redrafts", "lense", "lenes", "cantonal", "connatal", "curvidentate", "undetractive"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)