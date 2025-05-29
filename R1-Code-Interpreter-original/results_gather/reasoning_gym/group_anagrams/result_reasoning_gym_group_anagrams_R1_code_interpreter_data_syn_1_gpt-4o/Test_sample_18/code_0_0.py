from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["shouted", "southed", "chorions", "isochron", "desex", "dexes", "sexed", 
         "lilts", "tills", "still", "thein", "thine", "velicate", "celative", 
         "muratorian", "mortuarian", "fester", "freest", "tapery", "tepary", 
         "pratey", "petary"]

result = group_anagrams(words)
print(result)