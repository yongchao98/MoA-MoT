from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["ingresses", "singeress", "tumbles", "stumble", "groper", "porger", 
         "aesculin", "lunacies", "centrodorsal", "doroscentral", "dorsocentral", 
         "cloturing", "courtling", "pastimer", "primates", "ligroines", 
         "religions", "shames", "mashes"]

result = group_anagrams(words)
print(result)