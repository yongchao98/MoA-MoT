def group_anagrams(words):
    anagrams = {}
    for word in words:
        sorted_word = ''.join(sorted(word))
        if sorted_word in anagrams:
            anagrams[sorted_word].append(word)
        else:
            anagrams[sorted_word] = [word]
    return list(anagrams.values())

words = ["unstrict", "instruct", "ensient", "intense", "sennite", "sentine", 
         "markets", "estmark", "ailuro", "ourali", "sour", "ours", "rous", 
         "petroleur", "poulterer", "uglier", "guiler", "reguli", "gluier", 
         "ligure", "merosthenic", "microsthene", "surreal", "serrula"]

result = group_anagrams(words)
print(result)