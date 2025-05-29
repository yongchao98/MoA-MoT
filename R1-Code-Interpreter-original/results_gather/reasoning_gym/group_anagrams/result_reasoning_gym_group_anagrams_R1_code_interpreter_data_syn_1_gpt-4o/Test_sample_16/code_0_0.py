words = ["cp", "pc", "natalist", "atlantis", "horseless", "shoreless", "nonenergic", "nongeneric", "harateen", "theraean", "contoid", "odontic", "tampering", "premating", "gerundive", "ungrieved", "eslabon", "bonsela", "sonable", "soleness", "noseless"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)