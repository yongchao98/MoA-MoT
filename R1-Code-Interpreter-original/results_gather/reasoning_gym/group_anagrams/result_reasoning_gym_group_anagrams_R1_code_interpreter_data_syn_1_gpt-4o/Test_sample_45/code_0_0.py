from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["uncoif", "cuinfo", "browse", "bowers", "bowser", "clava", "caval", "shopbook", "bookshop", "silverer", "resilver", "sliverer", "revilers", "armor", "morra", "maror"]
result = group_anagrams(words)
print(result)