from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["taplash", "asphalt", "spathal", "melam", "lemma", "kirve", "kiver", "duali", "dulia", 
         "harpists", "tsarship", "starship", "wonnot", "wonton", "mowha", "mohwa", "caupo", 
         "acoup", "bucrnia", "rubican", "brucina"]

result = group_anagrams(words)
print(result)