words = ["yardland", "landyard", "peorian", "apeiron", "darns", "nards", "maternities", "semitertian", "seicento", "cenosite"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)