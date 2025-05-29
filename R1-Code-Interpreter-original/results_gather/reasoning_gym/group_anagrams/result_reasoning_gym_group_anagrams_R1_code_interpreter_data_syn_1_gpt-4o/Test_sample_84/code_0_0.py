words = ["peacoat", "opacate", "congaed", "decagon", "acephalan", "acalephan", "patte", "tapet", "longimanous", "longanimous", "haars", "raash", "haras", "asarh", "overturn", "turnover", "geomantic", "atmogenic", "fole", "floe", "fonder", "fronde"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)