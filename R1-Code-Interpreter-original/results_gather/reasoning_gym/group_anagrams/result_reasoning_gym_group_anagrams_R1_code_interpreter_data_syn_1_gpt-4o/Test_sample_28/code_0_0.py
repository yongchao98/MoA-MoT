words = ["discouple", "clupeoids", "infarction", "infraction", "draba", "barad", "chordal", "dorlach", "racers", "screar", "carers", "radiators", "sartoriad"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)