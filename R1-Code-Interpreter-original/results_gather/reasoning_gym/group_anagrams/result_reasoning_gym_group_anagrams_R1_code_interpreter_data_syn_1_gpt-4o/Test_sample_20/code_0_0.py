words = ["pialyn", "alypin", "plainy", "unsubtracted", "subtruncated", "prosecutory", "orycteropus", "legendist", "glistened", "fondu", "found", "shyster", "thyrses"]

from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = ''.join(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

result = group_anagrams(words)
print(result)