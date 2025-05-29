words = ["unmail", "lumina", "alumin", "alumni", "endothecia", "theodicean", "shipholder", "holdership", "mestees", "esteems", "semiprivate", "imperatives", "deduces", "seduced", "depeche", "cheeped"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)