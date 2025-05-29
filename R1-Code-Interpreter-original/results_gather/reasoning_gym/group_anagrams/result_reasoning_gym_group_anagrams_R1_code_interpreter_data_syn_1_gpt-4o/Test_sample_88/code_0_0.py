def group_anagrams(words):
    anagrams = {}
    for word in words:
        sorted_word = ''.join(sorted(word))
        if sorted_word in anagrams:
            anagrams[sorted_word].append(word)
        else:
            anagrams[sorted_word] = [word]
    return list(anagrams.values())

words = ["cerris", "criers", "ricers", "boarders", "reboards", "viatores", "travoise", "votaries", 
         "liberally", "libellary", "ovensman", "snemovna", "intertwin", "inwritten", "melonist", 
         "telonism", "xeronate", "exornate", "contadino", "anticodon", "cowardy", "cowyard"]

result = group_anagrams(words)
print(result)