def group_anagrams(words):
    anagrams = {}
    for word in words:
        sorted_word = ''.join(sorted(word))
        if sorted_word in anagrams:
            anagrams[sorted_word].append(word)
        else:
            anagrams[sorted_word] = [word]
    return list(anagrams.values())

words = ["drapers", "sparred", "ahmet", "mathe", "interplead", "peridental"]
result = group_anagrams(words)
print(result)