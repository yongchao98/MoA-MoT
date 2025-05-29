from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["prosely", "leprosy", "copras", "carpos", "siren", "reins", "rinse", "resin", "subarration", "saburration"]
result = group_anagrams(words)
print(result)