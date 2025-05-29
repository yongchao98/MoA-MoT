from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["merino", "inmore", "monier", "preacher", "repreach", "caurale", "arcuale", "diable", "deblai", "bailed", "baldie", "abdiel", "licorous", "couloirs", "sorgos", "grosso", "pharm", "prahm"]
result = group_anagrams(words)
print(result)