words = ["drisk", "dirks", "weatherbreak", "breakweather", "torsion", "isotron", "nitroso", "bikram", "imbark", "carisa", "sciara", "certifying", "rectifying", "tapir", "atrip", "pairt", "parti", "perchloric", "prechloric"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

print(group_anagrams(words))