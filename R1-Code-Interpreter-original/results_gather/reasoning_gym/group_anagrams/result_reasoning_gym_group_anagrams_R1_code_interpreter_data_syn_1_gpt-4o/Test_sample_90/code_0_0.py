words = ["caesuric", "curacies", "sailable", "isabella", "touristical", "staurolitic", "cooniest", "coonties", "contoise", "unpursuing", "unusurping", "somniate", "masonite", "outerness", "outreness", "salited", "details", "dilates", "distale", "praefects", "perfectas", "camoodie", "comoedia"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)