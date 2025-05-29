from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["incomputably", "uncompatibly", "wakeners", "rewakens", "airplot", "toprail", "portail", 
         "seedman", "demeans", "electromagnetic", "magnetoelectric", "monial", "nomial", 
         "almoin", "oilman", "fetichry", "fitchery"]

result = group_anagrams(words)
print(result)