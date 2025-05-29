from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["employers", "reemploys", "ciselure", "ciseleur", "crusilee", "poppel", "popple", "sepia", "paise", 
         "ropeman", "manrope", "proceeds", "procedes", "entrep", "repent", "terpen", "perten", "toodle", 
         "looted", "tooled", "evanid", "invade", "viande", "dread", "dared"]

result = group_anagrams(words)
print(result)