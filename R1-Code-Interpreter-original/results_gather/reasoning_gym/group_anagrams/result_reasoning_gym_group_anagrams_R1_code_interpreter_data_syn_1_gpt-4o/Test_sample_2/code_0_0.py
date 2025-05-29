from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["labiodental", "dentolabial", "intaglios", "ligations", "isolating", 
         "incorporated", "adrenotropic", "samiti", "timias", "ostracodan", 
         "cardosanto", "unadept", "untaped", "aeneid", "diaene", "aedine", 
         "trumperiness", "surprisement", "od", "do", "desmic", "medics"]

result = group_anagrams(words)
print(result)