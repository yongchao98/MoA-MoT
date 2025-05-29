from collections import defaultdict

def group_anagrams(words):
    anagrams = defaultdict(list)
    for word in words:
        sorted_word = tuple(sorted(word))
        anagrams[sorted_word].append(word)
    return list(anagrams.values())

words = ["dobe", "bode", "daleth", "lathed", "halted", "unredacted", "underacted", 
         "rebundle", "blendure", "latinos", "talions", "caried", "acider", "ericad", 
         "mugweed", "gumweed", "phaseouts", "taphouses"]

result = group_anagrams(words)
print(result)