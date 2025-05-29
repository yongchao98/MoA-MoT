from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["socotri", "octrois", "mus", "sum", "overglint", "revolting", "reluctant", "tralucent", 
         "calceus", "saccule", "moronities", "monitories", "rance", "nacre", "crena", "upwork", 
         "workup", "jarde", "jared"]

result = group_anagrams(words)
print(result)