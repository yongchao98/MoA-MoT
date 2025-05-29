from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["spignel", "spingel", "kra", "ark", "romance", "cremona", "ascebc", "ebcasc", "spoored", "prosode", "corrida", "ricardo", "brisk", "birks", "stract", "tracts", "tusseh", "tushes", "shutes", "calipers", "spiracle"]
result = group_anagrams(words)
print(result)