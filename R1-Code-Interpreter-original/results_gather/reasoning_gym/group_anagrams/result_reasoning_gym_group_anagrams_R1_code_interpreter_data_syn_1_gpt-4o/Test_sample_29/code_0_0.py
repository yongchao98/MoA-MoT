from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["brontides", "stonebird", "horotelic", "ochrolite", "lackers", "calkers", "slacker", "rarely", "rearly", "unhearsed", "unsheared", "criey", "ricey", "scomfit", "comfits", "pergamic", "crimpage", "unrip", "purin"]
result = group_anagrams(words)
print(result)