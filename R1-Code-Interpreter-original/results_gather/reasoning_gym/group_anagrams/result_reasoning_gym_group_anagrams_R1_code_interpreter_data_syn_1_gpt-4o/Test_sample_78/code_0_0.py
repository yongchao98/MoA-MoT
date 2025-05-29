words = ["preversing", "preserving", "reservist", "reversist", "sative", "itaves", "stevia", 
         "crinose", "sericon", "recoins", "unrebated", "eburnated", "underbeat", 
         "reflowers", "flowerers", "sermonical", "normalcies", "cisted", "edicts"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)